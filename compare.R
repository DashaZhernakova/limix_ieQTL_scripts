library(dplyr)

setwd("/Users/Dasha/work/UMCG/data/gender_differences/eQTLgen/test_limix/")
my_res <- read.delim("LLD_interactions_sex.my_R.txt", sep = "\t", as.is = T, check.names = F)
limix_res <- read.delim("LLD_interactions_sex.limix.txt", sep = "\t", as.is = T, check.names = F)

res_merged <-inner_join(my_res, limix_res, by = c("snp" = "snp_id", "gene" = "feature_id"))

pdf("../test_limix/limix_vs_my_R_comparison.pdf", height = 8, width = 6)
par(mfrow = c(2,1))
plot(res_merged$`p_geno:covariate`, res_merged$p_value, pch = 16, cex = 0.4, xlab = "Dasha's interaction P", ylab = "Marc Jan's interaction P")
plot(res_merged$`b_geno:covariate`, res_merged$beta,  pch = 16, cex = 0.4, xlab = "Dasha's interaction beta", ylab = "Marc Jan's interaction beta")
dev.off()
cor(abs(res_merged$`b_geno:covariate`), abs(res_merged$beta))
