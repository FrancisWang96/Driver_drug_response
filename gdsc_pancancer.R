library(tidyverse)
library(pheatmap)
library(wesanderson)
library(coin)

setwd("/Users/wangzeyuan/Library/Mobile Documents/com~apple~CloudDocs/gdsc")

gdsc_snv <- read.csv('mutations_summary_20230202.csv')
cgc1 <- read.csv('Census_allTue Nov 14 05_54_49 2023.csv')

# SKCM
drug_measures <- read.csv('PANCANCER_IC_Mon May 20 06_36_17 2024.csv')
cgc_tag <- 'melanoma'
drug_names <- c('dabrafenib', 'vemurafenib', 'dacarbazine', 'trametinib', 'cobimetinib', 'vindesine')

# snv auc ic50 
gdsc_snv_samples <- gdsc_snv$model_name%>%unique()
drug_measures_samples <- drug_measures$Cell.Line.Name%>%unique()
samples_intersect <- gdsc_snv_samples%>%intersect(drug_measures_samples)
gdsc_snv_sub <- gdsc_snv[gdsc_snv$model_name%in%samples_intersect,]
drug_measures_sub <- drug_measures[drug_measures$Cell.Line.Name%in%samples_intersect,]

gdsc_snv_sub_temp <- dplyr::select(gdsc_snv_sub, gene_symbol, model_name)

# 删掉重复
gdsc_snv_sub_temp$merged <- stringr::str_c(gdsc_snv_sub_temp$gene_symbol, gdsc_snv_sub_temp$model_name, sep = "_")
temp1 <- unique(gdsc_snv_sub_temp$merged)
temp2 <- stringr::str_split_fixed(temp1, "_", n=2)
temp2 <- dplyr::as_tibble(temp2)
names(temp2) <- names(gdsc_snv_sub_temp)

gene_symbol <- unique(temp2$gene_symbol)
model_name <- unique(temp2$model_name)

mutation_matrix <- matrix(0,length(model_name) , length(gene_symbol))
rownames(mutation_matrix) <- model_name
colnames(mutation_matrix) <- gene_symbol
for (i in 1:nrow(temp2)) {
  mutation_matrix[match(temp2$model_name[i],rownames(mutation_matrix)), match(temp2$gene_symbol[i],colnames(mutation_matrix))] <- 1
}

cgc_genes <- cgc1$Gene.Symbol[cgc1$Tumour.Types.Somatic.%>%str_detect(cgc_tag)]

if(!dir.exists('output')){
  dir.create('output')
}
setwd('output')

if(!dir.exists('SKCM')){
  dir.create('SKCM')
}
setwd('SKCM')

drug_names <- drug_names[drug_names%in%tolower(drug_measures_sub$Drug.Name)]
cgc_genes1 <- cgc_genes[cgc_genes%in%colnames(mutation_matrix)]

summary_auc <- data.frame(gene=rep(cgc_genes1,length(drug_names)),drug=rep(drug_names,each=length(cgc_genes1)), n_mutated=0,n_unmutated=0,p_ttest=2,p_mann=2,p_med=2)
summary_ic50 <- data.frame(gene=rep(cgc_genes1,length(drug_names)),drug=rep(drug_names,each=length(cgc_genes1)), n_mutated=0,n_unmutated=0,p_ttest=2,p_mann=2,p_med=2)


for (i in 1:length(cgc_genes1)) {
  for (drug in drug_names) {
    drug_measures_sub1 <- drug_measures_sub[tolower(drug_measures_sub$Drug.Name)==drug,]
    if(length(unique(drug_measures_sub1$Cell.Line.Name))==length(drug_measures_sub1$Cell.Line.Name)){
      mutation_matrix1 <- mutation_matrix[rownames(mutation_matrix)%in%drug_measures_sub1$Cell.Line.Name,]
      if(sum(mutation_matrix1[,match(cgc_genes1[i],colnames(mutation_matrix1))])>1){
        snv_auc <- data.frame(gene_symbol = mutation_matrix1[,match(cgc_genes1[i],colnames(mutation_matrix1))],
                              auc = drug_measures_sub1$AUC[match(drug_measures_sub1$Cell.Line.Name,rownames(mutation_matrix1))])
        snv_auc$gene_symbol[snv_auc$gene_symbol==1] <- str_c(cgc_genes1[i],'_mutated')
        snv_auc$gene_symbol[snv_auc$gene_symbol==0] <- str_c(cgc_genes1[i],'_unmutated')
        colnames(snv_auc)[1] <- 'Tag'
        snv_auc$Tag <- snv_auc$Tag%>%as.factor()
        snv_auc$auc <- scale(snv_auc$auc)
        
        summary_auc$n_mutated[summary_auc$gene==cgc_genes1[i]&summary_auc$drug==drug] <- sum(mutation_matrix1[,match(cgc_genes1[i],colnames(mutation_matrix1))])
        summary_auc$n_unmutated[summary_auc$gene==cgc_genes1[i]&summary_auc$drug==drug] <- nrow(mutation_matrix1)-sum(mutation_matrix1[,match(cgc_genes1[i],colnames(mutation_matrix1))])
        
        n_mutated <- sum(mutation_matrix1[,match(cgc_genes1[i],colnames(mutation_matrix1))])
        n_unmutated <- nrow(mutation_matrix1)-sum(mutation_matrix1[,match(cgc_genes1[i],colnames(mutation_matrix1))])
        
        auc_t_test_value <- t.test(snv_auc$auc ~ snv_auc$Tag)
        summary_auc$p_ttest[summary_auc$gene==cgc_genes1[i]&summary_auc$drug==drug] <- auc_t_test_value$p.value
        
        auc_w_test_value <- wilcox.test(snv_auc$auc ~ snv_auc$Tag)
        summary_auc$p_mann[summary_auc$gene==cgc_genes1[i]&summary_auc$drug==drug] <- auc_w_test_value$p.value
        
        summary_auc$p_med[summary_auc$gene==cgc_genes1[i]&summary_auc$drug==drug] <- pvalue(median_test(snv_auc$auc ~ snv_auc$Tag))
        
        
        # cat(t_test_value$p.value)
        # cat('\n')
        adjust=1
        snv_auc%>%ggplot(aes(y=auc))+
          geom_boxplot(aes(fill=Tag),color = 'black', alpha=0.4)+
          geom_text(aes(x=0,y=3+adjust,label=str_c('p.t = ',round(auc_t_test_value$p.value,3))))+
          geom_text(aes(x=0,y=2.6+adjust,label=str_c('p.w = ',round(auc_w_test_value$p.value,3))))+
          geom_text(aes(x=0,y=2.2+adjust,label=str_c('p.m = ',round(pvalue(median_test(snv_auc$auc ~ snv_auc$Tag)),3))))+
          geom_text(x=-0.2,y=2,label=str_c('n = ',n_mutated),color="#999999",size=5)+
          geom_text(x=0.2,y=2,label=str_c('n = ',n_unmutated),color="#E69F00",size=5)+
          scale_fill_manual(values=c("#999999", "#E69F00"))+
          theme_minimal()+
          labs(y ="AUC", x = "",title = str_c(cgc_genes1[i],' ',drug,' AUC'))
        ggsave(str_c('AUC_',cgc_genes1[i],'_',drug,'.pdf'),width=5,height = 4)
        
        snv_ic50 <- data.frame(gene_symbol = mutation_matrix1[,match(cgc_genes1[i],colnames(mutation_matrix1))],
                               ic50 = drug_measures_sub1$IC50[match(drug_measures_sub1$Cell.Line.Name,rownames(mutation_matrix1))])
        snv_ic50$gene_symbol[snv_ic50$gene_symbol==1] <- str_c(cgc_genes1[i],'_mutated')
        snv_ic50$gene_symbol[snv_ic50$gene_symbol==0] <- str_c(cgc_genes1[i],'_unmutated')
        colnames(snv_ic50)[1] <- 'Tag'
        snv_ic50$Tag <- snv_ic50$Tag%>%as.factor()
        snv_ic50$ic50 <- scale(snv_ic50$ic50)
        
        summary_ic50$n_mutated[summary_auc$gene==cgc_genes1[i]&summary_auc$drug==drug] <- sum(mutation_matrix1[,match(cgc_genes1[i],colnames(mutation_matrix1))])
        summary_ic50$n_unmutated[summary_auc$gene==cgc_genes1[i]&summary_auc$drug==drug] <- nrow(mutation_matrix1)-sum(mutation_matrix1[,match(cgc_genes1[i],colnames(mutation_matrix1))])
        
        t_test_value_ic50 <- t.test(snv_ic50$ic50 ~ snv_ic50$Tag)
        summary_ic50$p_ttest[summary_ic50$gene==cgc_genes1[i]&summary_ic50$drug==drug] <- t_test_value_ic50$p.value
        
        w_test_value_ic50 <- wilcox.test(snv_ic50$ic50 ~ snv_ic50$Tag)
        summary_ic50$p_mann[summary_ic50$gene==cgc_genes1[i]&summary_ic50$drug==drug] <- w_test_value_ic50$p.value
        
        summary_ic50$p_med[summary_ic50$gene==cgc_genes1[i]&summary_ic50$drug==drug] <- pvalue(median_test(snv_ic50$ic50 ~ snv_ic50$Tag))
        
        # cat(t_test_value$p.value)
        # cat('\n')
        adjust=1
        
        snv_ic50%>%ggplot(aes(y=ic50))+
          geom_boxplot(aes(fill=Tag),color = 'black', alpha=0.4)+
          geom_text(aes(x=0,y=3+adjust,label=str_c('p.t = ',round(t_test_value_ic50$p.value,3))))+
          geom_text(aes(x=0,y=2.6+adjust,label=str_c('p.w = ',round(w_test_value_ic50$p.value,3))))+
          geom_text(aes(x=0,y=2.2+adjust,label=str_c('p.m = ',round(pvalue(median_test(snv_ic50$ic50 ~ snv_ic50$Tag)),3))))+
          geom_text(x=-0.2,y=2,label=str_c('n = ',n_mutated),color="#999999",size=5)+
          geom_text(x=0.2,y=2,label=str_c('n = ',n_unmutated),color="#E69F00",size=5)+
          scale_fill_manual(values=c("#999999", "#E69F00"))+
          theme_minimal()+
          labs(y ="IC50", x = "",title = str_c(cgc_genes1[i],' ',drug,' IC50'))
        ggsave(str_c('IC50_',cgc_genes1[i],'_',drug,'.pdf'),width=5,height = 4)
        
      }
    }
  }
}
summary_auc <- summary_auc%>%arrange(p_ttest)
write_csv(summary_auc,'summary_auc.csv') 
summary_ic50 <- summary_ic50%>%arrange(p_ttest)
write_csv(summary_ic50,'summary_ic50.csv') 

setwd('..')
setwd('..')
