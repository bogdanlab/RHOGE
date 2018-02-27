
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
#' @return data.frame with estimates
get_cor <- function(t1, t2, regions = grch37.eur.loci, min_regions=10, nsamples=100) {
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
  out = list()
  for (idx in 1:nsamples) {
    sub_genes <- ind_genes %>% group_by(BLOCK) %>% sample_n(1) %>% pull(ID)
    sub_both <- both %>% filter(ID %in% sub_genes)
    out[idx] <- stats::cor(sub_both$Alpha.x, sub_both$Alpha.y)
  }
  vout <- unlist(out)
  m_rho_ge <- mean(vout)
  se_rho_ge <- sd(vout) / sqrt(nsamples)
  t_rho_ge <- m_rho_ge / se_rho_ge
  p <- 2 * pt(abs(t_rho_ge), mregions - 2, lower.tail=F)

  data.frame(RHOGE=m_rho_ge, SE=se_rho_ge, TSTAT=t_rho_ge, P=p)
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

  res1 <- get_trait(trait1, n1)
  res2 <- get_trait(trait2, n2)

  # Get nominally significant results both traits
  trt1 <- res1 %>% filter(TWAS.P < p)
  trt2 <- res2 %>% filter(TWAS.P < p)

  # Compute \rho_ge
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
#' @param min_genes, Minimum number of ascertained genes required for bi-directional regression
#' @param regions, data.frame-like containing approximate independent regions. Requires columns (Chr, Start, Stop). Default is estimated blocks in Europeans.
#' @export
rhoge.bd <- function(trait1, trait2, n1, n2, p1 = NA, p2 = NA, min_genes = 10, regions = grch37.eur.loci) {
  res1 <- get_trait(trait1, n1)
  res2 <- get_trait(trait2, n2)

  # Ascertain on genes specific to trait 1
  trt1 <- res1 %>% dplyr::filter(TWAS.P < p1)
  trt2 <- res2 %>% dplyr::filter(TWAS.P <= 1)

  res1 <- get_cor(trt1, trt2, regions = regions, min_regions = min_regions)
  rho12 <- res1$RHO
  p12 <- res1$P
  se12 <- res1$SE
  n12 <- res1$M

  # Ascertain on genes specific to trait 2
  trt1 <- res1 %>% dplyr::filter(TWAS.P <= 1)
  trt2 <- res2 %>% dplyr::filter(TWAS.P < p1)

  res2 <- get_cor(trt1, trt2, regions = regions, min_regions = min_regions)
  rho21 <- res2$RHO
  p21 <- res2$P
  se21 <- res2$SE
  n21 <- res2$M

  if (n21 > min_genes & n12 > min_genes) {
      # Compute if model means are different
      t <- (rho12 - rho21)/sqrt(se12^2 + se21^2)
      df <- (se12^2 + se21^2)^2/((se12^4/(n12 - 1)) + (se21^4/(n21 - 1)))
      Pt <- 2 * pt(abs(t), df, lower.tail = F)
  } else {
      df <- "NA"
      Pt <- "NA"
  }

  cat("(Trait1 -> Trait2) RhoGE SE P(t) M:", rho12, se12, p12, n12, "\n", file = output, append = T)
  cat("(Trait1 <- Trait2) RhoGE SE P(t) M:", rho21, se21, p21, n21, "\n", file = output, append = T)
  cat("P(Welch's t) ~df", Pt, df, "\n", file = output, append = T)
}
