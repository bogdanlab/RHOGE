library(argparser)
library(bootstrap)
library(plyr)
library(tidyr)


# Check input
check_args <- function(args) {
    return(TRUE)
}

# Get FUSION trait
get_trait <- function(path, NGWAS) {
    #FILE ID CHR P0 P1 HSQ BEST.GWAS.ID BEST.GWAS.Z EQTL.ID EQTL.R2 EQTL.Z EQTL.GWAS.Z NSNP NWGT MODEL MODELCV.R2 MODELCV.PV TWAS.Z TWAS.P PERM.PV PERM.N PERM.ANL_PV
    df <- read.table(path, h=T)

    # Filter on P-value
    df <- df[!is.na(df$TWAS.P) & !is.infinite(df$TWAS.Z) & !is.nan(df$TWAS.Z),]

    # Keep only...
    # FILE ID CHR P0 P1 HSQ TWAS.Z TWAS.P
    df <- df[,c("FILE", "ID", "CHR", "P0", "P1", "HSQ", "TWAS.Z", "TWAS.P")]

    # Remove SPLICING FOR NOW
    df <- df[!grepl("SPLICING", df$FILE),]

    # Get separate GENE table
    gene <- df[,c("ID", "CHR", "P0", "P1", "TWAS.P")]

    # Compute gene-expression effect, alpha
    df$Alpha <- df$TWAS.Z / sqrt(NGWAS * df$HSQ)
    df <- df[,c("FILE", "ID", "Alpha")]

    # Remove duplicates
    df$KEY <- paste0(df$FILE, df$ID)
    df <- df[!duplicated(df$KEY),]
    df$KEY <- NULL

    return (list(TWAS=df, GENES=gene))
}

get_table <- function(df) {
    ncol <- length(unique(df$FILE))

    cor_df <- spread(df, FILE, Alpha)
    cor_df$Alpha <- rowMeans(cor_df[,2:ncol + 1], na.rm=T)
    cor_df <- cor_df[,c("ID", "Alpha")]

    return (cor_df)
}

get_trait_genes <- function(genes, P) {
    sub <- genes[genes$TWAS.P < P,]
    sub <- sub[order(sub$TWAS.P),]
    sub <- sub[!duplicated(sub$ID),]
    return (sub)
}

# Compute RhoGE 
get_cor <- function(t1, t2, gene1, gene2, span) {

    ng1 <- nrow(gene1)
    ng2 <- nrow(gene2)
    if (is.null(ng1) | is.null(ng2) | ng1 == 0 | ng2 == 0) {
        return(list(RHO=NA, P=NA, SE=NA, M=0))
    }

    t1 <- get_table(t1[t1$ID %in% gene1$ID,])
    t2 <- get_table(t2[t2$ID %in% gene2$ID,])

    both <- merge(t1, t2, by="ID")

    # focus only on genes that have measurements in both traits
    both = both[!is.na(both$Alpha.x) & !is.na(both$Alpha.y),]

    # Gene sets are the same at this point, so use gene1 to grab info
    genes <- gene1[gene1$ID %in% both$ID,]

    # Grab the independent genes
    ind_genes <- get_independent(genes, span)
    M <- nrow(ind_genes)

    if (!is.null(M) & M > 0) {
        # Filter set to independent sites
        both = both[both$ID %in% ind_genes$ID,]

        M = nrow(both)
        # Compute correlation over independent genes
        RHO <- cor(both$Alpha.x, both$Alpha.y)

        # Jack-knife for standard error
        jck <- jackknife(1:nrow(both), function(x, xdata){cor(xdata[x,]$Alpha.x, xdata[x,]$Alpha.y, use="complete.obs")}, both)
        SE <- jck$jack.se

        # t ~ t(M - 2)
        P <- 2 * pt(abs(RHO / SE), M-2, lower.tail=F)
        return(list(RHO=RHO, P=P, SE=SE, M=M))
    } else {
        return(list(RHO=NA, P=NA, SE=NA, M=0))
    }
}

# Get spatially independent genes
get_independent <- function(genes, span) {
    if (is.null(genes) | nrow(genes) == 0) {
        return (c());
    }
    genes = genes[order(genes$CHR, genes$P0),]

    out = NULL
    for (chr in 1:22) {
        k = 1
        sub = genes[genes$CHR == chr,]
        if (dim(sub)[1] == 0) {
            next
        }
        if (is.null(out)) {
            out = as.data.frame(sub[1,])
        } else {
            out = rbind(out, sub[1,])
        }
        if (nrow(sub) == 1) {
            next
        }
        for (idx in 2:nrow(sub)) {
            if ( sub[idx,]$P1 >= (sub[k,]$P0 + span)) {
                out = rbind(out, sub[idx,])
                k = idx
            }
        }
    }
    return (out);
}

ap <- arg_parser("RhoGE: compute the genome-wide genetic correlation of two complex traits at the level of predicted (cis) expression")
ap <- add_argument(ap, "trait1", help="FUSION output for trait 1")
ap <- add_argument(ap, "trait2", help="FUSION output for trait 2")
ap <- add_argument(ap, "N1", type="integer", help="Sample size for trait 1")
ap <- add_argument(ap, "N2", type="integer", help="Sample size for trait 2")
ap <- add_argument(ap, "--P", type="double", default=0.05, help="Nominal significance threshold for initial RhoGE estimate.")
ap <- add_argument(ap, "--P1", type="double", default=0.0, help="Transcriptome-wide significance threshold for trait 1 ascertainment (Default is Bonferroni adjusted 0.05)")
ap <- add_argument(ap, "--P2", type="double", default=0.0, help="Transcriptome-wide significance threshold for trait 2 ascertainment (Default is Bonferroni adjusted 0.05)")
ap <- add_argument(ap, "--multi_tissue", short="-t", default=TRUE, flag=TRUE, help="Use the mean effect across tissues")
ap <- add_argument(ap, "--min_genes", short="-m", default=10, help="Required number of genes to perform bi-directional estimation")
ap <- add_argument(ap, "--region_size", short="-s", default=1e6, help="Size of genomic region when selecting independent genes")
ap <- add_argument(ap, "--output", short="-o", default="", help="Output to store results")

args <- parse_args(ap)

if (!check_args(args)) {
    quit(save="no", status=1)
}

res <- get_trait(args$trait1, args$N1)
trt1 <- res$TWAS
gene1 <- res$GENE
P1 <- ifelse (args$P1 != 0, args$P1, 0.05 / nrow(trt1))

if (is.null(trt1)) {
    quit(save="no", status=1)
}

res <- get_trait(args$trait2, args$N2)
trt2 <- res$TWAS
gene2 <- res$GENE
P2 <- ifelse (args$P2 != 0, args$P2, 0.05 / nrow(trt2))

if (is.null(trt2)) {
    quit(save="no", status=1)
}

# Get genes that are nominally significant for both traits
t1_genes <- get_trait_genes(gene1, args$P)
t2_genes <- get_trait_genes(gene2, args$P)

# Compute RHO_GE
res <- get_cor(trt1, trt2, t1_genes, t2_genes, args$region_size) 
if (!is.null(res)) {
    rho <- res$RHO
    se <- res$SE
    p <- res$P
    n <- res$M
} else {
    quit(save="no", status=1)
}

# Output before continuing...
cat("RhoGE SE P(t) M:", rho, se, p, n, "\n", file=args$output)

# Ascertain on genes specific to trait 1
t1_genes <- get_trait_genes(gene1, P1)
t2_genes <- get_trait_genes(gene2, 1)

res1 <- get_cor(trt1, trt2, t1_genes, t2_genes, args$region_size)
rho12 <- res1$RHO
p12 <- res1$P
se12 <- res1$SE
n12 <- res1$M

# Ascertain on genes specific to trait 2
t1_genes <- get_trait_genes(gene1, 1)
t2_genes <- get_trait_genes(gene2, P2)

res2 <- get_cor(trt1, trt2, t1_genes, t2_genes, args$region_size)
rho21 <- res2$RHO
p21 <- res2$P
se21 <- res2$SE
n21 <- res2$M

if (n21 > args$min_genes & n12 > args$min_genes) {
    # Compute if model means are different
    t <- (rho12 - rho21) / sqrt(se12^2 + se21^2)
    df <- (se12^2 + se21^2)^2 / ((se12^4 / (n12 - 1)) + (se21^4 / (n21 - 1)))
    Pt <- 2 * pt(abs(t), df, lower.tail=F)
} else {
    df <- "NA"
    Pt <- "NA"
}

cat("(Trait1 -> Trait2) RhoGE SE P(t) M:", rho12, se12, p12, n12, "\n", file=args$output, append=T)
cat("(Trait1 <- Trait2) RhoGE SE P(t) M:", rho21, se21, p21, n21, "\n", file=args$output, append=T)
cat("P(Welch's t) ~df", Pt, df, "\n", file=args$output, append=T)
