Sys.setenv(OPENGWAS_JWT="")
ieugwasr::get_opengwas_jwt()
ieugwasr::user()
library(data.table)
library(plyr)
library(dplyr)  
library(data.table)
library(TwoSampleMR)
library(readxl)
library(plinkbinr)
library(ClumpMR)

gene <- read_excel("gene.xlsx")
freq=fread("2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz")
a=gene$Ensembl
gene_name=gene$symbol


eqtl<-fread("2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded(2).txt.gz")
eqtl=subset(eqtl,GeneSymbol%in%gene_name)
gene_eqtl1=eqtl


head(freq)
freq1=freq[,c(1,9)]
gene_eqtl1=merge(gene_eqtl1,freq1,by="SNP")
colnames(gene_eqtl1)[16] <- "eaf"
gene_eqtl1$beta <- gene_eqtl1$Zscore / sqrt(2 * gene_eqtl1$eaf* (1 - gene_eqtl1$eaf) * (gene_eqtl1$NrSamples + gene_eqtl1$Zscore^2))
gene_eqtl1$SE <- abs(gene_eqtl1$beta/gene_eqtl1$Zscore)
gene_eqtl1$Zscore=NULL
gene_eqtl1$AlleleB_all <- NULL
data_clumped=list()
rm(gene)
b=gene_eqtl1$GeneSymbol
b=unique(b)

get_plink_exe()
for (i in 1:63) {
  data_not_clump=subset(gene_eqtl1,GeneSymbol==b[i])
  data_not_clump $id <- data_not_clump$GeneSymbol
  data_not_clump$rsid <- data_not_clump$SNP
  data_not_clump$pval <- data_not_clump$Pvalue
 data_clumped[[i]]=ieugwasr::ld_clump_local(data_not_clump,clump_r2=0.1,clump_kb=10000,clump_p=1,
                                            bfile="g1000_eur", 
                                            plink_bin="D:/R/R language/R-4.3.2/library/plinkbinr/bin/plink_Windows.exe")

}

data_exp=list()
for (i in 1:63) {
  data_clumped[[i]]=data.frame(data_clumped[[i]])
  data_exp[[i]]=data.frame(data_clumped[[i]]) %>% format_data(,type="exposure",
              phenotype_col = b[i],snp_col = "SNP",
              beta_col = "beta",se_col = "SE",
              pval_col = "Pvalue",
              samplesize_col = "NrSamples",
              effect_allele_col = "AssessedAllele",
              other_allele_col = "OtherAllele",
              chr_col = "SNPChr",
              pos_col = "SNPPos" )
}
results <- list()
resdata <- list()
res_hete <- list()
res_plei <- list()
res_presso <- list()
exp_data<- list()
out_come="Ulcerative colitis"


for (i in 1:63) {
  UC<-extract_outcome_data(
    snps = data_exp[[i]]$SNP,
    outcomes ='ebi-a-GCST004133',
    proxies = F,
    maf_threshold = 0.01,
    opengwas_jwt = ieugwasr::get_opengwas_jwt())
 
  dat <- harmonise_data(
    exposure_dat=data_exp[[i]],
    outcome_dat=UC,
    action= 2)
  dat <- subset(dat,mr_keep)
  res <- mr(dat,method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median"))
  OR<-generate_odds_ratios(res)
  OR
  results[[i]]=OR
  resdata[[i]] <- dat 
  if (nrow(dat[dat$mr_keep==T,])>1) {
    #sensitive analysis
    res1<- mr_heterogeneity(dat)
    res1$exposure=res$exposure[1]
    res1$outcome=out_come
    
    res2 <- mr_pleiotropy_test(dat)
    res2$exposure=res$exposure[1]
    res2$outcome=out_come
    
    res_hete[[i]] <- res1
    res_plei[[i]] <- res2
  }
  
  if (nrow(dat[dat$mr_keep==T,])>3) {
    set.seed(123)
    dd <- run_mr_presso(dat)
    res3 <- data.frame(exposure=res$exposure[1],outcome=out_come,pval=dd[[1]][["MR-PRESSO results"]][["Global Test"]][["Pvalue"]])
    res_presso[[i]] <- res3} 
  
}


results_allIV  <- do.call(rbind, results)
resdata_allIV <- do.call(rbind, resdata)
res_hete_allIV <- do.call(rbind, res_hete)
res_plei_allIV <- do.call(rbind, res_plei)
res_presso_allIV <- do.call(rbind, res_presso)

#format(round(x, 2), nsmall = 2)
results_allIV$estimate <- paste0(format(round(results_allIV$or, 2), nsmall = 2), "(", 
                                 format(round(results_allIV$or_lci95, 2), nsmall = 2),"-",
                                 format(round(results_allIV$or_uci95, 2), nsmall = 2), ")")




###################### output ###############
write.csv(results_allIV,file = "gene_results.csv")
write.csv(resdata_allIV,file = "gene_data.csv")
write.csv(res_hete_allIV,file = "gene_het.csv")
write.csv(res_plei_allIV,file = "gene_pleio.csv")
write.csv(res_presso_allIV,file = "gene_presso.csv")
save.image(file="final.RData")


  


