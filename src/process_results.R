
## ===============================================
##
## KIR gene PheWAS in FinnGen R11
## KIR-HLA interaction PheWAS in FinnGen R11
## Sum stat results plottting
##
## ===============================================

library(data.table)
library(tidyverse)
library(readxl)

library(ghibli)
library(wesanderson)
library(ggpubr)
library(scales)
library(patchwork)
library(ggbeeswarm)


## -----------------------------------------------
## Process summary stats
## -----------------------------------------------

# read sum stats
res <- fread('results/results_18Apr2024.tsv', data.table=F)

# pheno information
fg.pheno <- read_xlsx('data/FINNGEN_ENDPOINTS_DF11_Final_2022-10-05_public.xlsx') %>% 
  dplyr::select(c(TAGS, LEVEL, NAME, LONGNAME, HD_ICD_10, HD_ICD_9))
res <- left_join(res, fg.pheno, by=c('Pheno'='NAME'), keep = T)
# pheno group
res$Phenogroup <- str_split_fixed(res$Pheno, '_', 2)[, 1]
  
# rearrange cols
res <- res[, c(1, 24:25, 27:30, 2:4, 6, 8, 11:23)]

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
  out <- grepl(x, phenogroups$TAGS)
  if(all(out==F) | all(is.na(out))) NA else phenogroups$LONGNAME[which(out)[1]]
}) %>% unlist
res$phenogroup[is.na(res$phenogroup)] <- 'Other'

res <- res[sample(1:nrow(res), nrow(res), replace = F), ]
res <- arrange(res, phenogroup)
res$phenogroup <- factor(res$phenogroup, levels=res$phenogroup %>% unique)
res$Pheno <- factor(res$Pheno, levels=res$Phenoe %>% unique)


## -----------------------------------------------
## Plot significant assocs
## -----------------------------------------------

# KIR significance
res.signif <- res %>% filter(FDR_BY<0.05, TEST=='ADD')

# inlude sex
res.signif.sex <- map_dfr(1:nrow(res.signif), function(i) {
  filter(res, LONGNAME==res.signif$LONGNAME[i], 
         ID==res.signif$ID[i])
})

# re-format genotype names
res.signif.sex$ID <- gsub('C1_', 'C1-', res.signif.sex$ID) %>% gsub('KIR_', '       ', .)
# unite pheno and genotype names
res.signif.sex$ID <- str_pad(res.signif.sex$ID, 0, 'both')
res.signif.sex <- unite(res.signif.sex, PhenoID, LONGNAME, ID, sep='           ', remove=F)
# order by beta
res.signif.sex$PhenoID <- factor(res.signif.sex$PhenoID, 
                                 levels = res.signif.sex %>% filter(TEST=='ADD') %>% 
                                   group_by(PhenoID) %>% 
                                   summarise(PMean=min(Beta)) %>% 
                                   arrange(PMean) %>% .$PhenoID)
# significance variable
res.signif.sex$IsSignif <- ifelse(res.signif.sex$FDR_BY < 0.05, 'FDR <0.05', 'Not significant')

# forest plot
p.res.signif.sex <- ggplot(res.signif.sex, 
                           aes(Beta, PhenoID, xmax=CI95upper, xmin=CI95lower, color=TEST, shape=IsSignif)) +
  geom_pointrange(position=position_dodge(width = 0.7), fatten = 3.5) +
  scale_shape_manual(name = NULL, values = c(19, 1)) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = 'grey0', linetype = 'dashed') +
  coord_cartesian(clip = 'off', xlim = c(-1.3, 1.8), ylim = c(0.5, 13.3), expand=0) + 
  annotate('text', -2.72, 14.5, label = 'Phenotype') +
  annotate('text', -1.6, 14.5, label = 'Genotype') +
  scale_color_manual(values = c('red3', '#482878FF'),# c('salmon', hue_pal()(3)[c(3)]),  
                     labels = c('Genotype', 'Sex (male)'), name = NULL) +
  ggtitle('') +
  ylab('') +
  theme_minimal() +
  theme(axis.text.y=element_text(hjust = 1),
        panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank(),
        axis.line = element_line(linewidth = 0.3),
        legend.position='bottom')

p.res.signif.sex.n <- ggplot(res.signif.sex, aes(nCases/2, PhenoID)) +
  geom_bar(stat = 'identity', fill = 'grey', width = 0.6) +
  ggtitle('') +
  annotate('text', 1, 14, label='') +
  geom_vline(xintercept = 0, linewidth = 0.2, color = 'grey') +
  coord_cartesian(clip = 'off', ylim = c(0.5, 13.3), expand = 0) +
  ylab('') + xlab('# Cases') +
  geom_vline(xintercept = 10000, linewidth = 0.3, color = 'grey60', linetype = 'dashed') +
  geom_vline(xintercept = 20000, linewidth = 0.3, color = 'grey60', linetype = 'dashed') +
  geom_vline(xintercept = 30000, linewidth = 0.3, color = 'grey60', linetype = 'dashed') +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.line = element_line(linewidth = 0.3),
        legend.position = 'bottom')



## -----------------------------------------------
## plot manhattan
## -----------------------------------------------

pal <- c('#EABE94','#0B775E','#35274A','#F1401F','#899DA4','#C93312','#DC863B','#1F262EFF','#AE93BEFF',
         '#B4DAE4F1','#6FB382FF','#4D6D93FF','#4C413FFF','#132E41FF','#DE7862FF','#E75B64FF','#D8AF39FF',
         '#B6D1F241','#67B9E9FF','#DCCA2CFF','#E8C4A2FF','#5A6F80FF','#833437FF','#278B9AFF','#14191FFF')

p1 <- res %>% filter(TEST == 'ADD') %>% mutate(Signif = ifelse(FDR_BY<0.05, 'y', 'n')) %>% 
  ggplot(aes(-log10(P), phenogroup, color=phenogroup, alpha = Signif, shape = Signif)) +
  geom_jitter() +
  scale_alpha_manual(values = c(0.3, 1)) +
  scale_shape_manual(values = c(1, 19)) +
  scale_color_manual(values=pal) +
  geom_vline(xintercept = -log10(res %>% filter(TEST == 'ADD', FDR_BY<0.05) %>% .$P %>% max()), 
             linewidth = 0.3, color = 'grey0', linetype = 'dashed') +
  xlab(expression('-log'[10]*'('*italic(p)*')')) +
  ylab('Phenotype class') +
  coord_cartesian(clip='on', xlim = c(0, 4), expand = 0, ylim = c(0.5, 24.5)) +
  theme_minimal() +
  theme(axis.text.y=element_text(hjust=1),
        panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank(),
        axis.line = element_line(linewidth = 0.3),
        legend.position='none')



## -----------------------------------------------
## plot beta distr
## -----------------------------------------------

res$ID %>% unique
res2 <- res %>% filter(TEST == 'ADD') %>% group_by(ID, FDR_BY) %>% 
  summarise(M = mean(Beta), P = mean(-log10(FDR_BY)))
res2$ID <- gsub('C1_', 'C1-', res2$ID) %>% gsub('Bw4_', 'Bw4-', .) %>% gsub('KIR_', '       ', .)
res2$ID <- factor(res2$ID, 
                  levels = res2 %>% group_by(ID) %>% summarise(MM = mean(M)) %>% arrange(MM) %>% .$ID)

p2 <- res2 %>% ggplot(aes(M, ID, color = P)) +
  geom_boxplot(outliers = F, linewidth = 0.3) +
  geom_quasirandom(shape = 1, size = 0.8, alpha = 0.8) +
  scale_color_viridis_c(direction = -1, name = expression('-log'[10]*'('*italic(FDR)*')')) +
  coord_cartesian(clip = 'off') +
  ylab('Genotype') + xlab('Beta') +
  theme_minimal() +
  geom_vline(xintercept = 0, linewidth = 0.3, color = 'grey0', linetype = 'dashed') +
  theme(axis.text.y=element_text(hjust=1),
        panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank(),
        axis.line = element_line(linewidth = 0.3),
        legend.key.height = unit(0.5, 'cm'),
        legend.title.position = 'top',
        legend.position='top')


## -----------------------------------------------
## composite plot
## -----------------------------------------------

pdf('results/Figure1.pdf', width=11, height=9.5)
((p1 + p2 + plot_layout(widths = c(1, 0.68))) /
    (p.res.signif.sex + p.res.signif.sex.n + plot_layout(widths=c(1, 0.48)))
) + plot_layout(heights = c(1, 0.52)) +
  plot_annotation(tag_levels = list(c('a', 'b', 'c'))) & 
  theme(plot.tag = element_text(face = 'bold', size = 14))
dev.off()  

