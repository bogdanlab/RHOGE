
#' Internal parameter validation
#'
#' @param n size
#' @param name name of variable
validate_size <- function(n, name) {
  if (!is.numeric(n) || n < 0) {
    stop(paste0(name, " must be a positive integer value"))
  }
}

#' Internal pval validation
#' @param p double, nominal significance threshold for genome-wide RhoGE estimate (default=0.05)
#' @param name name of variable
validate_pval <- function(p, name) {
  if (!is.numeric(p) || !(p >= 0 && p <= 1)) {
    stop(paste0(name, " must be a number in [0, 1]"))
  }
}

#' Internal data.frame/table validation
#' @param trait1 data.frame-like, containing TWAS results for trait 1
#' @param trait2 data.frame-like, containing TWAS results for trait 2
#' @param regions, data.frame-like containing approximate independent regions. Requires columns (Chr, Start,
validate_tables <- function(trait1, trait2, regions = grch37.eur.loci) {
  check_df <- function(df, name) {
    if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) {
      stop(paste0(name, " must be a non-empty data.frame like object"))
    }
  }
  check_fusion_df <- function(df, name) {
    req_names <- c("FILE", "ID", "CHR", "P0", "P1", "HSQ", "TWAS.Z", "TWAS.P")
    cn <- colnames(df)
    rnn <- length(req_names)
    if (length(base::intersect(req_names, cn)) != rnn) {
      stop(paste0(name, " requires ", paste(req_names, sep=",", collapse = T), "for valid inference."))
    }
  }
  check_df(trait1, "trait1")
  check_fusion_df(trait1, "trait1")
  check_df(trait2, "trait2")
  check_fusion_df(trait2, "trait2")

  check_df(regions, "regions")
  if (!("CHR" %in% colnames(regions) && "START" %in% colnames(regions) && "STOP" %in% colnames(regions))) {
    stop("regions must be non-empty data.frame-like with CHR, START, and STOP variables")
  }
}

#' Get FUSION table
#' @param twas_tbl data.frame-like, Table/Dataframe of TWAS results from FUSION
#' @param ngwas integer, Sample size for GWAS
#' @importFrom dplyr filter select mutate
#' @return Annotated data.frame suitable for RHOGE inference
get_trait <- function(twas_tbl, ngwas) {
  # FUSION-style header FILE ID CHR P0 P1 HSQ BEST.GWAS.ID BEST.GWAS.Z EQTL.ID EQTL.R2 EQTL.Z EQTL.GWAS.Z NSNP NWGT MODEL MODELCV.R2 MODELCV.PV TWAS.Z TWAS.P
  twas_tbl <- twas_tbl %>% filter(!is.na(TWAS.P) & !is.infinite(TWAS.Z) & !is.nan(TWAS.Z))

  # Keep only...  FILE ID CHR P0 P1 HSQ TWAS.Z TWAS.P
  twas_tbl <- twas_tbl %>% select(FILE, ID, CHR, P0, P1, HSQ, TWAS.Z, TWAS.P)

  # Compute gene-expression effect, alpha
  twas_tbl %>% mutate(Alpha = TWAS.Z/sqrt(ngwas * HSQ), KEY = paste0(FILE, ID)) %>%
    filter(!duplicated(KEY)) %>% select(FILE, ID, CHR, P0, P1, Alpha, TWAS.P)
}

#' Compute \eqn{\rho_ge}
#' @param t1 data.frame-like, TWAS results for trait1
#' @param t2 data.frame-like, TWAS results for trait2
#' @param regions data.frame-like, Approximate independent regions. Requires columns (Chr, Start, Stop). Default is estimated blocks in Europeans.
#' @param min_regions integer, Minimum number of independent regions for inference
#' @param nsamples integer, Number of samples to perform for bootstrapping
#' @importFrom dplyr group_by summarize n_distinct pull filter sample_n
#' @importFrom bootstrap jackknife
#' @return data.frame with estimates
get_cor <- function(t1, t2, regions = grch37.eur.loci, min_regions=10, nsamples=10) {
  # mean is likely over-confident since predicted expression levels will be moderately correlated across tissues
  # TODO: consider other meta-analysis methods in future
  t1 <- t1 %>% group_by(ID) %>% summarize(CHR=CHR[1], P0=median(P0), P1=median(P1), Alpha=mean(Alpha))
  t2 <- t2 %>% group_by(ID) %>% summarize(Alpha=mean(Alpha))

  # focus only on genes that have measurements in both traits
  both <- merge(t1, t2, by = "ID") %>% filter(!is.na(Alpha.x) & !is.na(Alpha.y))

  # Grab the independent genes
  ind_genes <- get_independent(both, regions)
  mregions <- n_distinct(ind_genes$BLOCK)

  if (mregions < min_regions) {
    stop(paste0("Only estimated ", mregions, " independent regions. Minimum number of regions required for reliable estimation is 10."))
  }

  # bootstrap genes from each each region
  # this turns out to be very stable and produces biased standard errors
  # to remedy this just do a small amount of bootstraps and jackknife the sd each time
  cor.out = list()
  sd.out = list()
  for (idx in 1:nsamples) {
    sub_genes <- ind_genes %>% group_by(BLOCK) %>% sample_n(1) %>% pull(ID)
    sub_both <- both %>% filter(ID %in% sub_genes)
    cor.out[idx] <- stats::cor(sub_both$Alpha.x, sub_both$Alpha.y)

    jkres <- bootstrap::jackknife(1:length(sub_genes), function(x) stats::cor(sub_both[x,]$Alpha.x, sub_both[x,]$Alpha.y))
    sd.out[idx] <- jkres$jack.se
  }
  vout <- unlist(cor.out)
  vsdout <- unlist(sd.out)

  m_rho_ge <- mean(vout)
  se_rho_ge <- mean(vsdout)

  t_rho_ge <- m_rho_ge / se_rho_ge
  p <- 2 * pt(abs(t_rho_ge), mregions - 2, lower.tail=F)

  data.frame(RHOGE=m_rho_ge, SE=se_rho_ge, TSTAT=t_rho_ge, DF=mregions - 2, P=p)
}

#' Get spatially independent genes
#'
#' @param genes data.frame-like, gene information for a TWAS
#' @param regions data.frame-like, independent regions contained in Chr, Start, Stop variables
#' @importFrom dplyr mutate filter
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom IRanges IRanges
#' @return ind_genes data.frame-like, data.frame annotated with independent blocks for each gene
get_independent <- function(genes, regions) {
  if (is.null(genes) || nrow(genes) == 0 ||
      !("CHR" %in% colnames(genes) && "P0" %in% colnames(genes) && "P1" %in% colnames(genes))) {
    stop("genes must be non-empty data.frame-like with CHR, P0, and P1 variables")
  }
  if (is.null(regions) || nrow(regions) == 0 ||
      !("CHR" %in% colnames(regions) && "START" %in% colnames(regions) && "STOP" %in% colnames(regions))) {
    stop("regions must be non-empty data.frame-like with CHR, START, and STOP variables")
  }

  # I really don't like depending on so much machinery for such a simple overlap check...
  g_ranges <- GenomicRanges::GRanges(seqnames = genes$CHR, ranges= IRanges::IRanges(start = genes$P0, end = genes$P1, names = genes$ID))
  r_ranges <- GenomicRanges::GRanges(seqnames = regions$CHR, ranges= IRanges::IRanges(start = regions$START, end = regions$STOP))
  ovlaps <- GenomicRanges::findOverlaps(g_ranges, r_ranges, select="arbitrary")

  # Sample a random gene from the region
  genes %>% dplyr::mutate(BLOCK=ovlaps) %>% dplyr::filter(!is.na(BLOCK))
}

#' Calculates genome-wide genetic correlation between two complex traits using TWAS summary statistics.
#'
#' @param trait1 data.frame-like, TWAS results for trait 1
#' @param trait2 data.frame-like, TWAS results for trait 2
#' @param n1 integer, sample-size for trait 1 GWAS
#' @param n2 integer, sample-size for trait 2 GWAS
#' @param p double, nominal significance threshold for genome-wide RhoGE estimate (default=0.05)
#' @param regions data.frame-like, Approximate independent regions. Requires columns (Chr, Start, Stop). Default is estimated blocks in Europeans.
#' @importFrom dplyr filter
#' @export
rhoge.gw <- function(trait1, trait2, n1, n2, p = 0.05, regions = grch37.eur.loci) {

  # perform error-checking here
  validate_tables(trait1, trait2, regions)
  validate_size(n1, "n1")
  validate_size(n2, "n2")
  validate_pval(p, "p")

  res1 <- get_trait(trait1, n1)
  res2 <- get_trait(trait2, n2)

  # Get nominally significant results both traits
  trt1 <- res1 %>% filter(TWAS.P < p)
  trt2 <- res2 %>% filter(TWAS.P < p)

  # Compute \eqn{\rho_ge}
  get_cor(trt1, trt2, regions)
}

#' TBD
#'
#' @param trait1 data.frame-like, containing TWAS results for trait 1
#' @param trait2 data.frame-like, containing TWAS results for trait 2
#' @param n1 integer, sample-size for trait 1 GWAS
#' @param n2 integer, sample-size for trait 2 GWAS
#' @param p1 double, Transcriptome-wide significance threshold for trait 1 ascertainment. Default is # Bonferroni adjusted 0.05
#' @param p2 double, Transcriptome-wide significance threshold for trait 2 ascertainment. Default is # Bonferroni adjusted 0.05
#' @param min_regions, Minimum number of ascertained regions required for bi-directional regression
#' @param regions, data.frame-like containing approximate independent regions. Requires columns (Chr, Start, Stop). Default is estimated blocks in Europeans.
#' @importFrom dplyr filter rename mutate bind_rows
#' @export
rhoge.bd <- function(trait1, trait2, n1, n2, p1 = NA, p2 = NA, min_regions = 10, regions = grch37.eur.loci) {
  # perform error-checking here
  validate_tables(trait1, trait2, regions)
  validate_size(n1, "n1")
  validate_size(n2, "n2")
  validate_size(min_regions, "min_regions")

  res1 <- get_trait(trait1, n1)
  res2 <- get_trait(trait2, n2)

  if (is.na(p1)) {
    p1 <- 0.05 / nrow(trait1)
  } else {
    validate_pval(p1, "p1")
  }

  if (is.na(p2)) {
    p2 <- 0.05 / nrow(trait2)
  } else {
    validate_pval(p2, "p2")
  }

  # Ascertain on genes specific to trait 1
  trt1 <- res1 %>% dplyr::filter(TWAS.P < p1)
  trt2 <- res2 %>% dplyr::filter(TWAS.P <= 1)

  rhoge1 <- get_cor(trt1, trt2, regions = regions, min_regions = min_regions)
  rho12 <- rhoge1$RHOGE
  p12 <- rhoge1$P
  se12 <- rhoge1$SE
  n12 <- rhoge1$DF + 2

  # Ascertain on genes specific to trait 2
  trt1 <- res1 %>% dplyr::filter(TWAS.P <= 1)
  trt2 <- res2 %>% dplyr::filter(TWAS.P < p2)

  rhoge2 <- get_cor(trt1, trt2, regions = regions, min_regions = min_regions)
  rho21 <- rhoge2$RHOGE
  p21 <- rhoge2$P
  se21 <- rhoge2$SE
  n21 <- rhoge2$DF + 2

  f1 <- rhoge1 %>% dplyr::rename(ESTIMATE = RHOGE) %>% dplyr::mutate(TEST = "Trait1 -> Trait2")
  f2 <- rhoge2 %>% dplyr::rename(ESTIMATE = RHOGE) %>% dplyr::mutate(TEST = "Trait2 -> Trait1")

  if (n21 > min_regions & n12 > min_regions) {
    t <- (rho12 - rho21)/sqrt(se12^2 + se21^2)
    df <- (se12^2 + se21^2)^2/((se12^4/(n12 - 1)) + (se21^4/(n21 - 1)))
    Pt <- 2 * pt(abs(t), df, lower.tail = F)
    diff <- data.frame(ESTIMATE=NA, SE=NA, TEST="DIFF", TSTAT=t, DF=df, P=Pt)

  } else {
    warning("Insufficient observations to compute difference test.")
    diff <- data.frame(ESTIMATE=NA, SE=NA, TEST="DIFF", TSTAT=NA, DF=NA, P=NA)
  }

  bind_rows(list(f1, f2, diff))
}
