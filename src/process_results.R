
## ===============================================
##
## KIR gene PheWAS in FinnGen R11
## KIR-HLA interaction PheWAS in FinnGen R11
##
## ===============================================

library(data.table)
library(tidyverse)
library(readxl)

library(ghibli)
library(wesanderson)
library(ggpubr)



## -----------------------------------------------
## Process summary stats
## -----------------------------------------------

# read sum stats
res <- fread('results/results_06May2023.tsv', data.table=F)
# re-orient HLA-KIR genotype results

# pheno information
fg.pheno <- read_xlsx('data/FINNGEN_ENDPOINTS_DF11_Final_2022-10-05_public.xlsx') %>% 
  dplyr::select(c(TAGS, LEVEL, NAME, LONGNAME, HD_ICD_10, HD_ICD_9))
res <- left_join(res, fg.pheno, by=c('Pheno'='NAME'))
# pheno group
res$Phenogroup <- str_split_fixed(res$Pheno, '_', 2)[, 1]
  
# rearrange cols
res <- res[, c(21, 1, 19, 24, 22:23, 2:18)]
colnames(res)[1:2] <- c('Phenotype', 'Phenotype short')
res$REF <- ifelse(res$REF=='P', 'Present', 'Absent')
res$ALT <- ifelse(res$ALT=='P', 'Present', 'Absent')
res$A1 <- ifelse(res$A1=='P', 'Present', 'Absent')

fwrite(res, 'results/summary_stats.tsv', sep='\t')

## -----------------------------------------------
## High-level phenotype grouping
## -----------------------------------------------

phenogroups <- filter(fg.pheno, LEVEL==1)
phenogroups$TAGS <- gsub('#', '', phenogroups$TAGS, fixed=T)
phenogroups$TAGS <- gsub('ICDMAIN,', '', phenogroups$TAGS, fixed=T)
phenogroups$TAGS <- gsub(',AVO|,OPHTAL|,GASTRO_CM|,CARDIO,NEURO_CM|,GASTRO', '', phenogroups$TAGS)
phenogroups <- phenogroups[, c('TAGS', 'LONGNAME')] %>% data.frame()
phenogroups <- phenogroups[-c(8, 11, 27, 28), ]
phenogroups$LONGNAME[8] <- 'In situ and benign neoplasms'

res$phenogroup <- map(res$TAGS %>% gsub('#', '', .) %>% gsub(',', '|', .), function(x) {
  # x <- "K11,DENTAL"
  # print(x)
  out <- grepl(x, phenogroups$TAGS)
  if(all(out==F) | all(is.na(out))) NA else phenogroups$LONGNAME[which(out)[1]]
}) %>% unlist
res$phenogroup[is.na(res$phenogroup)] <- 'Other'

res <- res[sample(1:nrow(res), nrow(res), replace = F), ]
res <- arrange(res, phenogroup)
res$phenogroup <- factor(res$phenogroup, levels=res$phenogroup %>% unique)
res$Phenotype <- factor(res$Phenotype, levels=res$Phenotype %>% unique)


## -----------------------------------------------
## plot
## -----------------------------------------------


# pal <- c(ghibli_palettes$SpiritedMedium, ghibli_palettes$PonyoMedium, 
#          ghibli_palettes$LaputaMedium, ghibli_palettes$YesterdayMedium)
# pal <- sample(pal, 25, replace = F)
# # pal.save <- pal
# pal2 <- pal.save
# pal2[1:7] <- c(wes_palettes$Rushmore, wes_palettes$Royal1)[c(2:7,9)]
# pal2

pal <- c('#EABE94','#0B775E','#35274A','#F2300F','#899DA4','#C93312','#DC863B','#1F262EFF','#AE93BEFF',
         '#B4DAE5FF','#6FB382FF','#4D6D93FF','#4C413FFF','#132E41FF','#DE7862FF','#E75B64FF','#D8AF39FF',
         '#B7D9F2FF','#67B9E9FF','#DCCA2CFF','#E8C4A2FF','#5A6F80FF','#833437FF','#278B9AFF','#14191FFF')

p1 <- ggplot(res, aes(Phenotype, mlog10P, color=phenogroup, fill=phenogroup)) +
  geom_jitter(size=1.25, shape=24, alpha=0.75, data=filter(res, Z_STAT>0)) +
  geom_jitter(size=1.25, shape=25, alpha=0.75, data=filter(res, Z_STAT<0)) +
  # geom_jitter(size=1.25, shape=21, alpha=0.7) +
  scale_color_manual(values=pal) +
  scale_fill_manual(values=pal) +
  ylab(expression('-log'[10]*'('*italic(p)*')')) +
  xlab('Phenotype class') +
  # scale_y_continuous(expand = expansion(mult=c(5.5, 5.5))) +
  coord_cartesian(clip='off', ylim = c(0, 5), expand=F) +
  facet_wrap(phenogroup ~., scales = 'free_x', nrow = 1, strip.position = 'bottom') +
  # geom_hline(yintercept = -0.058, color='white') +
  # theme_pubclean() +
  theme_minimal() +
  # geom_hline(yintercept = 0, color = 'black', linewidth = 0.4) +
  annotate(geom = 'line', 0.4, y=0) +
  theme(legend.position = 'none',
        axis.line.y = element_line(color = 'black', linewidth = 0.4),
        axis.line.x = element_line(color = 'white', linewidth = 0.4),
        axis.text.x = element_blank(),
        axis.title.x = element_text(margin = margin(-20, 0, 0, 0)),
        panel.spacing = unit(.1, 'lines'),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(angle = 25, hjust=1.0, vjust=1.0, size = 6.5),
        axis.ticks.x = element_blank(),
        panel.background = element_blank(),
        plot.margin = margin(5, 5, 10, 125, unit='pt'))

jpeg('results/Fig1.jpg', width = 9, height = 4.5, res = 600, units = 'in')
p1
dev.off()


## -----------------------------------------------
## pheno - ICD10 table
## -----------------------------------------------

res$phenogroup %>% unique %>% map_dfr(function(x) {
  out <- filter(res, phenogroup==x)$HD_ICD_10 %>% na.omit %>% unique %>% paste(collapse=', ')
  data.frame(Pheno_class=x, ICD10=out) %>% return
}) %>% fwrite(., 'results/phenotype_class_ICD10.tsv', sep='\t')




