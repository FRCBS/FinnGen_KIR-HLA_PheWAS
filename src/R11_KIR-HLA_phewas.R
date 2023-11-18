
## ===============================================
##
## KIR gene PheWAS in FinnGen R11
## KIR-HLA interaction PheWAS in FinnGen R11
##
## ===============================================


library(data.table)
library(tidyverse)


## -----------------------------------------------
## kinship
## -----------------------------------------------

kin <- fread('/finngen/library-red/finngen_R11/kinship_1.0/data/finngen_R11.kin0', data.table=F)
kin <- filter(kin, InfType=='FS' | InfType=='PO' | InfType=='Dup/MZ')
kin$ID2 # 1st degree relatives to be removed from data
kin <- kin[, c('ID2', 'FID2')]
write(kin[, 1], 'data/phenotypes/plink_remove_ID', ncolumns=1)


## -----------------------------------------------
## phenotypes
## -----------------------------------------------

# extract fenotype file
system(paste0("cp /finngen/library-red/finngen_R11/phenotype_1.0/data/finngen_R11_endpoint_1.0.txt.gz ", 
              "data/phenotypes/finngen_R11_endpoint_1.0.txt.gz"))
system(paste0("gunzip data/phenotypes/finngen_R11_endpoint_1.0.txt.gz"))

# list phenos
system("head -n 1 data/phenotypes/finngen_R11_endpoint_1.0.txt > data/phenotypes/pheno.header")
pheno.header <- fread('data/phenotypes/pheno.header', data.table=F) %>% colnames
# phenotype selection
pheno.index <- which(!grepl('_NEVT|_AGE|_YEAR', pheno.header))
pheno.index <- pheno.index[-c(2:3)]
pheno.index <- pheno.index %>% paste(collapse=',')

# select columns from pheno file
system(paste0("cut -f", pheno.index, " data/phenotypes/finngen_R11_endpoint_1.0.txt > data/phenotypes/R11_cleaned_phenos.tsv"))

# remove original large pheno file
#system("rm data/phenotypes/finngen_R11_endpoint_1.0.txt")

# read cleaned phenos
phenos <- fread('data/phenotypes/R11_cleaned_phenos.tsv', data.table=F)
phenos %>% dim
# number of cases per pheno
phenos.n <- phenos[, -1] %>% colSums()
phenos.n %>% hist(100)
# filter out too rare and too common
phenos.keep <- which(phenos.n > 300 & phenos.n < 2e5)+1
phenos.keep %>% length
phenos <- phenos[, c(1, phenos.keep)]
phenos %>% dim()
colnames(phenos)[1] <- 'IID'
phenos <- data.frame(FID=phenos$IID, phenos)
# write plink pheno file
fwrite(phenos, 'data/phenotypes/R11_plink_phenos.tsv', sep='\t')
# test pheno
tmp <- fread('data/phenotypes/plink_test_pheno', data.table=F)
tmp <- tmp[sample(1:nrow(tmp), nrow(tmp), replace=F), ]
fwrite(tmp, 'data/phenotypes/plink_test_pheno2', sep='\t')

# n cases per pheno
pheno.ncases <- filter(phenos, !(IID %in% kin$ID2))[, -c(1:2)] %>% colSums()
pheno.ncases <- data.frame(Pheno=names(pheno.ncases), 
                           nCases=pheno.ncases)
fwrite(pheno.ncases, 'data/phenotypes/R11_plink_phenos_ncases.tsv', sep='\t')


# HLA ligands
# 2DL2, 2DL1, 3DL1, 3DL2
ligand.groups <- list(
  C1='C\\*01\\:02|C\\*03\\:02|C\\*03\\:03|C\\*03\\:04|C\\*07\\:01|C\\*07\\:02|C\\*07\\:04|C\\*08\\:02|C\\*12\\:02|C\\*12\\:03|C\\*14\\:02|C\\*16\\:01',
  C2='C\\*02\\:02|C\\*04\\:01|C\\*04\\:06|C\\*05\\:01|C\\*06\\:02|C\\*15\\:02|C\\*15\\:05|C\\*16\\:02|C\\*17\\:01|C\\*17\\:03',
  Bw4_80I='A\\*23\\:01|A\\*24\\:02|A\\*24\\:07|A\\*25\\:01|A\\*32\\:01|B\\*27\\:02|B\\*38\\:01|B\\*49\\:01|B\\*51\\:01|B\\*57\\:01',
  A3_A11='A\\*03\\:01|A\\*11\\:01'
)




## -----------------------------------------------
## covariates
## -----------------------------------------------

covars <- fread('data/phenotypes/R11_covar.txt', data.table=F)
covars %>% colnames
# covars <- dplyr::select(covars, -c(batch, n_var, chip, AGE_AT_DEATH_OR_END_OF_FOLLOWUP2, BL_YEAR,
#                                    BL_AGE, regionofbirth, COHORT, 145:154, 157:158, 160:166, 177:188)) # 114:141, 14:113
covars <- dplyr::select(covars, c(FID, IID, AGE_AT_DEATH_OR_END_OF_FOLLOWUP, SEX, IS_AFFY, regionofbirthname, 167:176))
# covars <- data.frame(FID=covars$IID, covars)
covars <- covars %>% na.omit
fwrite(covars, 'data/phenotypes/R11_covar_plink.tsv', sep='\t')



## -----------------------------------------------
## KIR genotypes
## -----------------------------------------------

kir.dat <- map2(list.files('results/KIR_imputation', 'imputed', full.names=T), 
               list.files('results/KIR_imputation', 'imputed', full.names=F) %>% gsub('imputed_|.tsv', '', .),
               function(x, y) {
                 out <- fread(x, data.table=F)[, -2] # this selects the PP for gene presence
                 colnames(out)[2] <- y
                 out %>% return
                }) %>% Reduce(left_join, .)

# func for plink conversion
toPlinkDos <- function(x, out.path.name) { 
  # input data: col 1 = sample ID, cols 2..n = geno dosage
  x[, -1] <- 2-x[, -1] 
  tmp <- x %>% t
  tmp.sam <- tmp[1, ]
  tmp.sam <- paste(tmp.sam, tmp.sam)
  tmp <- tmp[-1, ] %>% data.frame
  tmp <- data.frame(SNP=rownames(tmp),
                    A1='A',
                    A2='P',
                    tmp)
  colnames(tmp)[-c(1:3)] <- tmp.sam
  tmp.map <- data.frame(str_split_fixed(tmp.sam, ' ', 2), 0, 0, 0, 0)
  fwrite(tmp, 'tmp/tmp.dosage', sep='\t')
  fwrite(tmp.map, 'tmp/tmp.dosage.map', sep='\t', col.names=F)
  
  system(paste0("plink2 --import-dosage tmp/tmp.dosage single-chr=19 format=1 ",
                "--fam tmp/tmp.dosage.map --make-pgen --out ", out.path.name))
}

# write plink dosage file
toPlinkDos(kir.dat, 'data/genotypes/R11_KIR')
# check
system("plink2 --pfile data/genotypes/R11_KIR --recode A --out data/genotypes/R11_KIR_dos")

# read plink dos, compare with kir.dat to check orientation
# the plink dos file records gene presence
tmp <- fread('data/genotypes/R11_KIR_dos.raw')


## -----------------------------------------------
## HLA genotypes
## -----------------------------------------------

# HLA BFEN file
hla.bgen   <- '/finngen/library-red/finngen_R11/hla_1.0/R11_HLA.bgen'
hla.sample <- '/finngen/library-red/finngen_R11/hla_1.0/R11_HLA.bgen.sample'

# convert to plink dosage format
system(paste0("plink2 --bgen ", hla.bgen, " ref-first --sample ", hla.sample, 
              " --sort-vars --make-pgen --out data/genotypes/R11_HLA"))
system(paste0("plink2 --pfile data/genotypes/R11_HLA --recode A --out data/genotypes/R11_HLA_dos"))

# HLA plink dosage file records allele absence
hla <- fread('data/genotypes/R11_HLA_dos.raw', data.table=F)[, -c(3:6)]
# convert to allele presence
hla[, -c(1:2)] <- 2-hla[, -c(1:2)] # present dosage
colnames(hla)[-c(1:2)] <- gsub('_<absent>', '', colnames(hla)[-c(1:2)])

# class I
hla <- hla[, c(1:2, which(grepl('^A|^B|^C', colnames(hla)[-c(1:2)])))]

# hla ligand group dosages
hla.c1  <- hla[, c(1:2, which(grepl(ligand.groups[[1]], colnames(hla))))]
hla.c2  <- hla[, c(1:2, which(grepl(ligand.groups[[2]], colnames(hla))))]
hla.bw4 <- hla[, c(1:2, which(grepl(ligand.groups[[3]], colnames(hla))))]
hla.a3  <- hla[, c(1:2, which(grepl(ligand.groups[[4]], colnames(hla))))]

# calculate group dosage as maximum of all HLAs in the group
hla.c1  <- data.frame(hla.c1[,  c(1:2)], C1=hla.c1[,  -c(1:2)] %>% apply(., 1, max))
hla.c2  <- data.frame(hla.c2[,  c(1:2)], C2=hla.c2[,  -c(1:2)] %>% apply(., 1, max))
hla.bw4 <- data.frame(hla.bw4[, c(1:2)], Bw4=hla.bw4[,-c(1:2)] %>% apply(., 1, max))
hla.a3  <- data.frame(hla.a3[,  c(1:2)], A3=hla.a3[,  -c(1:2)] %>% apply(., 1, max))


## -----------------------------------------------
## HLA-KIR genotypes
## -----------------------------------------------

# 2DL2
hla.c1.2DL2 <- inner_join(hla.c1, kir.dat[, c('sample.id', 'KIR_2DL2')], by=c('IID'='sample.id'))
hla.c1.2DL2$C1_2DL2 <- hla.c1.2DL2$C1*hla.c1.2DL2$KIR_2DL2

# 2DL1
hla.c2.2DL1 <- inner_join(hla.c2, kir.dat[, c('sample.id', 'KIR_2DL1')], by=c('IID'='sample.id'))
hla.c2.2DL1$C2_2DL1 <- hla.c2.2DL1$C2*hla.c2.2DL1$KIR_2DL1

# 3DL1
hla.bw4.3DL1 <- inner_join(hla.bw4, kir.dat[, c('sample.id', 'KIR_3DL1')], by=c('IID'='sample.id'))
hla.bw4.3DL1$Bw4_3DL1 <- hla.bw4.3DL1$Bw4*hla.bw4.3DL1$KIR_3DL1

# 3DL2 is not imputed
hla.a3.3DL2 <- inner_join(hla.a3, kir.dat[, c('sample.id', 'KIR_3DL2')], by=c('IID'='sample.id'))
hla.a3.3DL2$A3_3DL2 <- hla.a3.3DL2$A3*hla.a3.3DL2$KIR_3DL2

# func for plink conversion
toPlinkDos2 <- function(x, out.path.name) { 
  # input data: col 1 = sample ID, cols 2..n = geno dosage
  # x <- hla.c1.2DL2[, c('IID', 'C1_2DL2')] %>% head
  # x <- hla.c1.2DL2[, c('IID', 'C1', 'C1_2DL2')] %>% head
  x[, -1] <- 2 - x[, -1] # switch orientation
  x <- data.frame(x, dummy=rep(0, nrow(x)))
  tmp <- x %>% t
  tmp.sam <- tmp[1, ]
  tmp.sam <- paste(tmp.sam, tmp.sam)
  tmp <- tmp[-1, ] %>% data.frame
  tmp <- data.frame(SNP=rownames(tmp),
                    A1='A',
                    A2='P',
                    tmp)
  colnames(tmp)[-c(1:3)] <- tmp.sam
  tmp <- tmp[-nrow(tmp), ]
  tmp.map <- data.frame(str_split_fixed(tmp.sam, ' ', 2), 0, 0, 0, 0)
  fwrite(tmp, 'tmp/tmp.dosage', sep='\t')
  fwrite(tmp.map, 'tmp/tmp.dosage.map', sep='\t', col.names=F)
  
  system(paste0("plink2 --import-dosage tmp/tmp.dosage single-chr=19 format=1 ",
                "--fam tmp/tmp.dosage.map --make-pgen --out ", out.path.name))
}


# write plink dosage files
toPlinkDos2(hla.c1.2DL2[, c('IID', 'C1_2DL2')], 'data/genotypes/R11_C1_2DL2')
toPlinkDos2(hla.c2.2DL1[, c('IID', 'C2_2DL1')], 'data/genotypes/R11_C2_2DL1')
toPlinkDos2(hla.bw4.3DL1[, c('IID', 'Bw4_3DL1')], 'data/genotypes/R11_Bw4_3DL1')

# check
system(paste0("plink2 --pfile data/genotypes/R11_C1_2DL2 --recode A --out tmp/tmp_dos"))
tmp <- fread('tmp/tmp_dos.raw')
hla.c1.2DL2[, c('IID', 'C1_2DL2')] %>% head


## -----------------------------------------------
## HLA-KIR covariate files
## -----------------------------------------------

covars.c1.2DL2  <- inner_join(covars, hla.c1.2DL2[,  c('IID', 'C1',  'KIR_2DL2')])
fwrite(covars.c1.2DL2, 'data/phenotypes/R11_covar_plink_c1.2DL2.tsv', sep='\t')

covars.c2.2DL1  <- inner_join(covars, hla.c2.2DL1[,  c('IID', 'C2',  'KIR_2DL1')])
fwrite(covars.c2.2DL1, 'data/phenotypes/R11_covar_plink_c2.2DL1.tsv', sep='\t')

covars.bw4.3DL1 <- inner_join(covars, hla.bw4.3DL1[, c('IID', 'Bw4', 'KIR_3DL1')])
fwrite(covars.bw4.3DL1, 'data/phenotypes/R11_covar_plink_bw4.3DL1.tsv', sep='\t')


## -----------------------------------------------
## Association
## -----------------------------------------------

# test
system(paste0("plink2 --pfile data/genotypes/R11_KIR --remove-fam data/phenotypes/plink_remove_ID.tsv ", 
              "--glm hide-covar --ci 0.95  --covar-variance-standardize ",
              "--covar data/phenotypes/R11_covar_plink.tsv ",
              "--pheno data/phenotypes/plink_test_pheno2 --1 ", 
              "--out tmp/test_KIR_phewas2"))

# phewas KIR
system(paste0("plink2 --pfile data/genotypes/R11_KIR --remove-fam data/phenotypes/plink_remove_ID.tsv ", 
              "--glm hide-covar --ci 0.95  --covar-variance-standardize ",
              "--covar data/phenotypes/R11_covar_plink.tsv ",
              "--pheno data/phenotypes/R11_plink_phenos.tsv --1 ", 
              "--out results/KIR_phewas/KIR_phewas"))

# phewas HLA C1 KIR 2DL2
system(paste0("plink2 --pfile data/genotypes/R11_C1_2DL2 --remove-fam data/phenotypes/plink_remove_ID.tsv ", 
              "--glm hide-covar --ci 0.95  --covar-variance-standardize ",
              "--covar data/phenotypes/R11_covar_plink_c1.2DL2.tsv ",
              "--pheno data/phenotypes/R11_plink_phenos.tsv --1 ", 
              "--out results/KIR_phewas/C1_2DL2_phewas"))

# phewas HLA C2 KIR 2DL1
system(paste0("plink2 --pfile data/genotypes/R11_C2_2DL1 --remove-fam data/phenotypes/plink_remove_ID.tsv ", 
              "--glm hide-covar --ci 0.95  --covar-variance-standardize ",
              "--covar data/phenotypes/R11_covar_plink_c2.2DL1.tsv ",
              "--pheno data/phenotypes/R11_plink_phenos.tsv --1 ", 
              "--out results/KIR_phewas/C2_2DL1_phewas"))

# phewas HLA Bw4 KIR 3DL1
system(paste0("plink2 --pfile data/genotypes/R11_Bw4_3DL1 --remove-fam data/phenotypes/plink_remove_ID.tsv ", 
              "--glm hide-covar --ci 0.95  --covar-variance-standardize ",
              "--covar data/phenotypes/R11_covar_plink_bw4.3DL1.tsv ",
              "--pheno data/phenotypes/R11_plink_phenos.tsv --1 ", 
              "--out results/KIR_phewas/Bw4_3DL1_phewas"))


## -----------------------------------------------
## Results
## -----------------------------------------------

# read reseult summaries
res <- map2(list.files('results/KIR_phewas', 'glm', full.names=T),
            list.files('results/KIR_phewas', 'glm', full.names=F) %>% 
              gsub('KIR_phewas.|C1_2DL2_phewas.|Bw4_3DL1_phewas.|C2_2DL1_phewas|.glm.logistic.hybrid', '', .), 
            function(x, y) {
              out <- fread(x, data.table=F)
              data.frame(Pheno=y, out) %>% return
            }) %>% do.call(rbind, .) %>% arrange(P)
res <- res %>% na.omit
res$mlog10P <- -log10(res$P)
res$FDR_BY <- p.adjust(res$P, method='BY')

# include case n
pheno.ncases <- fread('data/phenotypes/R11_plink_phenos_ncases.tsv', data.table=F)
res <- left_join(res, pheno.ncases, by=c('Pheno'='Pheno'))
res <- res[, c(1, 4:10, 20, 11:19)] # rearrange cols

# check
res$nCases %>% min() # 219
res %>% head(30)
fwrite(res, 'results/KIR_phewas/results_06May2023.tsv', sep='\t')



