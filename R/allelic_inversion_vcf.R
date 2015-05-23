################################################# ALLELIC INVERSION VCF #################################################
# #get the info from the vcf file in the Terminal run this within the working directory
grep -v "^#" batch_1.vcf | cut -f 2,4,5,8 | perl -pe 's/;/\t/g' | perl -pe 's/,/\t/g' | perl -pe 's/NS=//g' | perl -pe 's/AF=//g' > vcf.ref.allele

# back in R
vcf.ref.allele <- read.table("vcf.ref.allele",header=F)
colnames(vcf.ref.allele) <- c("POS", "REF", "ALT", "N_VCF", "FREQ_REF", "GLOBAL_MAF")

sumstats.inversion.vcf <- vcf.ref.allele %>%
  mutate(POS = POS-1) %>%
  select (POS, REF, N_VCF, GLOBAL_MAF) %>%
  merge (
    (arrange (sumstats, POS)), by="POS"
    ) %>%
  mutate(
    INVERSION = ifelse(as.character(ALLELE_P) == as.character(REF),"no","inversion"),
    FREQ_REF = ifelse(INVERSION == "inversion", FREQ_ALLELE_Q, FREQ_ALLELE_P),
    FREQ_ALT = ifelse(INVERSION == "inversion", FREQ_ALLELE_P, FREQ_ALLELE_Q)
    )%>%
  arrange(LOCUS, POS, POP_ID)

# Get the number of inversion...
number.inversion <- sumstats.inversion.vcf%>%
  filter (INVERSION == "inversion") %>%
  select (POS, POP_ID)

n_distinct(number.inversion$POS)
number.inversion.pop <- data.frame(xtabs(data=number.inversion,~POP_ID))

# reorder sumstats
sumstats.inversion.vcf<-sumstats.inversion.vcf[c("BATCH","LOCUS","CHROM","POS","COL","POP_ID","N", "N_VCF", "ALLELE_P", "REF", "INVERSION","ALLELE_Q","FREQ_ALLELE_P", "FREQ_REF", "FREQ_ALLELE_Q", "FREQ_ALT", "GLOBAL_MAF","HET_O","HOM_O","HET_E","HOM_E","PI","SMOOTHED_PI","SMOOTHED_PI_P_VALUE","FIS","SMOOTHED_FIS","SMOOTHED_FIS_P_VALUE","PRIVATE")]

#save the result:
write.table(sumstats.inversion.vcf,"sumstats.inversion.vcf",sep="\t",row.names=F,col.names=T,quote=F)

# to re-import:
sumstats.inversion.vcf <-read.table ("sumstats.inversion.vcf", header=T, dec=".", sep="\t")

# clean your desk...
rm(number.inversion.pop, number.inversion, vcf.ref.allele)
