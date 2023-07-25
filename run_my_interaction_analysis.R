args <- commandArgs(trailingOnly = TRUE)
library(dplyr)
library(patchwork)
library(ieugwasr)
library(ggplot2)

setwd("/Users/Dasha/work/UMCG/data/gender_differences/eQTLgen/v3/")
cohort="LLD"
covariate = "sex"
expr_fname <- paste0("data/", cohort, ".gene_read_counts_BIOS_and_LLD_passQC.TMM.ProbesCentered.SamplesZTransformed.txt.gz")
geno_fname <- paste0("data/", cohort, ".eqtl_genotypes.filtered.dosages.txt.gz")
gte_fname <- paste0("data/", cohort, "_samples.txt")
covar_fname <- paste0("cell_counts/tmp/", cohort, ".TPM.nnls.LM22+sex.txt")
eqtls_fname <- "data/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.top_hits.txt"
out_fname <- "../test_limix/LLD_interactions_sex.my_R.txt"

#covariate = "sex"
is_binary = ifelse(covariate == "sex", T, F)

# Read in the data
expr_noadj <- as.data.frame(t(read.delim(gzfile(expr_fname), check.names = F, header = T, row.names = 1)))
geno <- as.data.frame(t(read.delim(gzfile(geno_fname), check.names = F, header = T, row.names = 1)))
gte <- read.delim(gte_fname, check.names = F, header = F)
sex <- read.delim(covar_fname, check.names = F, header = T, row.names = 1)
sex <- na.omit(sex)

expr_noadj <- apply(expr_noadj, 2, function(x) qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x))))
contin_covars <- apply(sex, 2, function(x) length(unique(x)) > 3)
sex_int <- cbind(sex[,!contin_covars], apply(sex[,contin_covars], 2, function(x) qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x))) ))
colnames(sex_int) <- c(colnames(sex)[!contin_covars], colnames(sex)[contin_covars])

# rename expression sample ids to genotype sample ids
geno_m <- geno[match(gte[,2], row.names(geno), nomatch = 0),]
gte_m <- gte[match(row.names(geno_m), gte[,2], nomatch = 0),]
row.names(geno_m) <- gte_m[,3]
cat("Coverted genotype sample ids to expression sample ids.\nN samples before: ", nrow(geno), "\nN samples after: ", nrow(geno_m), "\n")

# get ids for samples shared between expr_noadjession,  genotypes and sex
ids <- intersect(intersect(row.names(expr_noadj), row.names(geno_m)), row.names(sex_int))
cat ("Number of overlapping samples between genotype and expr_noadjession data: ", length(ids), "\n")
expr_noadj2 <- expr_noadj[ids,]
geno2 <- geno_m[ids, ]
sex2 <- sex_int[ids,]
sex2 <- subset(sex2, select=-c(LDC))
colnames(sex2) <- gsub("gender_F1M2", "sex", colnames(sex2))

# Read eQTLs to test:
eqtls_full <- read.delim(eqtls_fname, check.names = F, header = T, as.is = T)
eqtls <- eqtls_full[eqtls_full$SNP %in% colnames(geno2),]
cat("Read ", nrow(eqtls), " eQTLs")

covars <- c("sex","age", "coravgexp","B cells naive", "Monocytes", "NK cells resting", "Neutrophils", "T cells CD8", "T cells regulatory (Tregs)")

# Run analysis
res <- data.frame(matrix(nrow = nrow(eqtls), ncol = 9))
colnames(res) <- c("gene", "snp", "b_geno", "b_covariate", "b_geno:covariate", "p_geno", "p_covariate", "p_geno:covariate", "se_geno:covariate")
cnt <- 1
for (e in 1:nrow(eqtls)){
  gene <- eqtls[e, "Gene"]
  snp <- eqtls[e, "SNP"]
  #cat(gene, snp,"\n")
  if (gene %in% colnames(expr_noadj2) & snp %in% colnames(geno2)){
    m <- as.data.frame(cbind(expr_noadj2[, gene], geno2[, snp], sex2[,covars]))
    colnames(m) <- c("gene", "dosage", covars)
    colnames(m) <- gsub("gender_F1M2", "sex", colnames(m))
    m$sex <- as.factor(m$sex)
    colnames(m) <- gsub(paste0('\\b', covariate, '\\b'), "covariate", colnames(m))
    row.names(m) <- ids
    genotype <- as.factor(round(m$dosage))
    
    if (min(table(genotype)) > 10){ # skip cases when # alt homo is less than 10
      
      lm_fit <- lm(gene ~ . + covariate*dosage , data = m)
      #summary(lm_fit)
      #af_sex <- wilcox.test(dosage ~ sex, data = m)$p.value
      coef <- summary(lm_fit)$coefficients
      row.names(coef) <- gsub("covariate2", "covariate", row.names(coef))
      
      res[cnt,] <- c(gene, snp, coef["dosage",1], coef["covariate",1], coef["dosage:covariate",1],  coef["dosage",4], coef["covariate",4], coef["dosage:covariate",4], coef["dosage:covariate",2])
      cnt <- cnt + 1
    }
  }
  if (cnt %% 10000 == 0){
    cat("processed ", cnt, " eQTLs\n")
  }
}

res[,3:(ncol(res)-1)] = apply(res[,3:(ncol(res)-1)], 2, function(x) as.numeric(as.character(x)))

res_srt <- res[order(res[,"p_geno:covariate"]),]
res_srt[,"p_geno:covariate_BH"] <- p.adjust(res_srt[,"p_geno:covariate"], method = "BH")
write.table(res_srt, file = out_fname, sep = "\t", quote = F, col.names = NA)


