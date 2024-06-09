library(tidyverse)
library(pheatmap)
library(wesanderson)
library(RColorBrewer)
library(viridis)
library(splines)


setwd("/Users/wangzeyuan/Library/Mobile Documents/com~apple~CloudDocs/ccle/23")

cell_line_info23 <- read.csv('Model.csv')

setwd("/Users/wangzeyuan/Library/Mobile Documents/com~apple~CloudDocs/ccle")

cell_line_info <- read.csv('secondary-screen-cell-line-info.csv')

cell_line_info <- cell_line_info%>%subset(!is.na(depmap_id))
cell_line_info$Model <- cell_line_info23$DepmapModelType[match(cell_line_info$depmap_id,cell_line_info23$ModelID)]

ccle_mutation <- read.csv('CCLE_mutations.csv')
ccle_cnv <- read.csv('CCLE_gene_cn.csv')
cgc <- read.csv('Census_allTue Nov 14 01_26_11 2023.csv')
cgc1 <- read.csv('Census_allTue Nov 14 05_54_49 2023.csv')
drug_response <- read.csv('secondary-screen-dose-response-curve-parameters.csv')
drug_response <- drug_response[!is.na(drug_response$indication),]

disease_list <- c('BLCA','GB','MEL','PAAD','LUAD','COAD')
cgc_list <- c('bladder','glioblastoma','melanoma','pancreatic','NSCLC','colon')
disease_list_drug <- c('bladder cancer','glioblastoma','melanoma','pancreatic cancer','non-small cell lung cancer','colorectal')

temp <- dplyr::select(ccle_mutation, Hugo_Symbol, DepMap_ID)
# 删掉重复
temp$merged <- stringr::str_c(temp$Hugo_Symbol, temp$DepMap_ID, sep = "_")
temp1 <- unique(temp$merged)
temp2 <- stringr::str_split_fixed(temp1, "_", n=2)
temp2 <- dplyr::as_tibble(temp2)
names(temp2) <- names(temp)
Hugo_Symbol <- unique(temp2$Hugo_Symbol)
Hugo_Symbol_table <- dplyr::tibble(Hugo_Symbol = Hugo_Symbol, flag = 1:length(Hugo_Symbol))

DepMap_ID <- unique(temp2$DepMap_ID)
DepMap_ID_table <- dplyr::tibble(DepMap_ID = DepMap_ID, flag = 1:length(DepMap_ID))

mutation_matrix <- matrix(0,length(DepMap_ID) , length(Hugo_Symbol))
row.names(mutation_matrix) <- DepMap_ID
colnames(mutation_matrix) <- Hugo_Symbol

temp2$r <- DepMap_ID_table$flag[match(temp2$DepMap_ID,
                                      DepMap_ID_table$DepMap_ID)]
temp2$c <- Hugo_Symbol_table$flag[match(temp2$Hugo_Symbol,
                                        Hugo_Symbol_table$Hugo_Symbol)]

for (i in 1:nrow(temp2)) {
  if(i %% floor(nrow(temp2)/10) == 0){
    cat(stringr::str_c(i %/% floor(nrow(temp2)/10),'0% '))
  }
  mutation_matrix[temp2$r[i],temp2$c[i]] <- 1
}

if(!dir.exists('single_nsclc')){
  dir.create('single_nsclc')
}
setwd('single_nsclc')

disease=5

if(!dir.exists(disease_list[disease])){
  dir.create(disease_list[disease])
}
setwd(disease_list[disease])

drug_response_breast <- drug_response[drug_response$indication%>%str_detect(disease_list_drug[disease]),]
cell_id <- cell_line_info$depmap_id[cell_line_info$Model==disease_list[disease]]
drug_response_breast <- drug_response_breast[drug_response_breast$depmap_id%in%cell_id,]
# flag <- drug_response_breast$depmap_id%>%match(ccle_expression$X)
# drug_response_breast <- drug_response_breast[!is.na(flag),]
drug_names <- drug_response_breast$name%>%unique()

drug_response_breast$auc_sd <- drug_response_breast$auc%>%scale()

if(!dir.exists('drugauc_distribution')){
  dir.create('drugauc_distribution')
}
setwd('drugauc_distribution')

drug_response_breast%>%ggplot(aes(x=auc_sd,..density..))+
  geom_histogram(colour="black", fill="white")+
  geom_density(alpha=0.2, fill = "coral3")+
  theme_minimal()+
  labs(title=str_c('Standardized General AUC Distribution on ',disease_list[disease]),x='AUC',y='Density')
ggsave('sd_auc_dis.pdf',width=5,height = 3)

for (i in drug_names) {
  drug_response_breast_sub <- drug_response_breast[drug_response_breast$name==i,]
  output_file <- str_c('sd_',i,'.pdf')
  drug_response_breast_sub%>%ggplot(aes(x=auc_sd,..density..))+
    geom_histogram(colour="black", fill="white")+
    geom_density(alpha=0.2, fill = "coral")+
    theme_minimal()+
    labs(title=str_c('Standardized AUC Distribution of ',i),x='AUC',y='Density')
  ggsave(output_file,width=5,height = 3)
}
setwd('..')

if(!dir.exists('drugec50_distribution')){
  dir.create('drugec50_distribution')
}
setwd('drugec50_distribution')

drug_response_breast1 <- drug_response_breast
drug_response_breast1$lgec50 <- drug_response_breast1$ec50%>%log()
drug_response_breast1$lgec50[drug_response_breast1$lgec50>20] <- 25
drug_response_breast1$lgec50 <- drug_response_breast1$lgec50%>%scale()
drug_response_breast$lgec50 <- drug_response_breast1$lgec50

drug_response_breast1%>%ggplot(aes(x=lgec50))+
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=0.5, fill = "aquamarine3")+
  theme_minimal()+
  labs(title=str_c('Standardized General log(EC50) Distributionon on ',disease_list[disease]),x='log(EC50)',y='Density')
ggsave('sd_ec50_dis.pdf',width=5,height = 3)

for (i in drug_names) {
  drug_response_breast_sub <- drug_response_breast[drug_response_breast$name==i,]
  drug_response_breast_sub$lgec50 <- drug_response_breast_sub$ec50%>%log()
  drug_response_breast_sub$lgec50[drug_response_breast_sub$lgec50>20] <- 25
  drug_response_breast_sub$lgec50 <- drug_response_breast_sub$lgec50%>%scale()
  output_file <- str_c('sd_',i,'.pdf')
  drug_response_breast_sub%>%ggplot(aes(x=lgec50))+
    geom_histogram(aes(y=..density..), colour="black", fill="white")+
    geom_density(alpha=0.5, fill = "aquamarine")+
    theme_minimal()+
    labs(title=str_c('Standardized log(EC50) Distribution of ',i),x='log(EC50)',y='Density')
  ggsave(output_file,width=5,height = 3)
}
setwd('..')

# drug_response_breast1%>%ggplot(aes(x=ec50))+
#   geom_histogram(aes(y=..density..), colour="black", fill="white")+
#   geom_density(alpha=0.2, fill = "#E69F00")+
#   theme_minimal()+
#   labs(y='Density')

# flag <- drug_response_breast$depmap_id%>%match(ccle_expression$X)
# ccle_expression <- ccle_expression[flag,]
# ccle_expression$drug <- drug_response_breast$name
# ccle_expression$ec50 <- drug_response_breast$ec50
# ccle_expression$auc <- drug_response_breast$auc
# ccle_expression$ic50 <- drug_response_breast$ic50
# 
# colnames(ccle_expression)

flag <- rownames(mutation_matrix)%in%drug_response_breast$depmap_id
mutation_matrix_breast <- mutation_matrix[flag,]

flag <- cgc1$Tumour.Types.Somatic.%>%str_detect(cgc_list[disease])
cgc_breast <- cgc1[flag,]
cgc_flag <- cgc_breast$Gene.Symbol%>%match(colnames(mutation_matrix_breast))
cgc_breast <- cgc_breast[!is.na(cgc_flag),]
cgc_flag <- cgc_flag[!is.na(cgc_flag)]
patient_cgc <- c()
patient_ncgc <- c()

for (i in 1:nrow(mutation_matrix_breast)) {
  if(sum(mutation_matrix_breast[i,cgc_flag])>0){
    patient_cgc <- patient_cgc%>%append(rownames(mutation_matrix_breast)[i])
  }else{
    patient_ncgc <- patient_ncgc%>%append(rownames(mutation_matrix_breast)[i])
  }
}

flag1 <- drug_response_breast$depmap_id%in%patient_cgc
drug_response_breast_cgc <- drug_response_breast[flag1,]
flag2 <- drug_response_breast$depmap_id%in%patient_ncgc
drug_response_breast_ncgc <- drug_response_breast[flag2,]
#
drug_response_breast_all <- drug_response_breast_cgc%>%rbind(drug_response_breast_ncgc)
drug_response_breast_all$tag <- c(rep('Driver_mutated',nrow(drug_response_breast_cgc)),rep('Driver_unmutated',nrow(drug_response_breast_ncgc)))

mutation_p <- data.frame(value = drug_response_breast_all$auc_sd, tag = drug_response_breast_all$tag, name = drug_response_breast_all$name)
mutaion_pvalues <- data.frame(drug=drug_names,p_t=1,p_w=1,p_m=1)
for (i in 1:length(drug_names)) {
  mutation_p1 <- mutation_p%>%subset(name==drug_names[i])
  mutation_p1$tag <- mutation_p1$tag%>%as.factor()
  t_test_value <- t.test(mutation_p1$value ~ mutation_p1$tag)
  mutaion_pvalues$p_t[i] <- t_test_value$p.value
  w_test_value <- wilcox.test(mutation_p1$value ~ mutation_p1$tag)
  mutaion_pvalues$p_w[i] <- w_test_value$p.value
  mutaion_pvalues$p_m[i] <- pvalue(median_test(mutation_p1$value ~ mutation_p1$tag))
}

dataframeforplot <- data.frame(value = drug_response_breast_all$auc_sd, tag = drug_response_breast_all$tag, name = drug_response_breast_all$name)


dataframeforplot1 <- dataframeforplot%>%distinct(name,tag)
dataframeforplot1$p_t <- round(mutaion_pvalues$p_t,3)
dataframeforplot1$p_w <- round(mutaion_pvalues$p_w,3)
dataframeforplot1$p_m <- round(mutaion_pvalues$p_m,3)
dataframeforplot1$sig.threshold <- as.factor(ifelse(dataframeforplot1$p_t<0.01,'p<0.01',ifelse(dataframeforplot1$p_t<0.05,'p<0.05',ifelse(dataframeforplot1$p_t<0.1,'p<0.1','p>=0.1'))))
dataframeforplot1$w.sig.threshold <- as.factor(ifelse(dataframeforplot1$p_w<0.01,'p<0.01',ifelse(dataframeforplot1$p_w<0.05,'p<0.05',ifelse(dataframeforplot1$p_w<0.1,'p<0.1','p>=0.1'))))
dataframeforplot1$m.sig.threshold <- as.factor(ifelse(dataframeforplot1$p_m<0.01,'p<0.01',ifelse(dataframeforplot1$p_m<0.05,'p<0.05',ifelse(dataframeforplot1$p_m<0.1,'p<0.1','p>=0.1'))))

dataframeforplot1$sig.threshold <- factor(dataframeforplot1$sig.threshold,levels =c('p<0.01','p<0.05','p<0.1','p>=0.1'))
dataframeforplot1$w.sig.threshold <- factor(dataframeforplot1$w.sig.threshold,levels =c('p<0.01','p<0.05','p<0.1','p>=0.1'))
dataframeforplot1$m.sig.threshold <- factor(dataframeforplot1$m.sig.threshold,levels =c('p<0.01','p<0.05','p<0.1','p>=0.1'))
adjust=-1.5

dataframeforplot%>%ggplot(aes(x = name,y=value))+
  geom_boxplot(aes(fill = tag),color = 'black', alpha=0.4)+
  stat_summary(fun = mean,geom = "point", shape=5, size=1.5,aes(group=tag), position = position_dodge(0.8))+
  geom_text(data=dataframeforplot1,aes(y=4.8+adjust,label=str_c('p.t=',p_t),color=sig.threshold),size=3)+
  geom_text(data=dataframeforplot1,aes(y=4.6+adjust,label=str_c('p.w=',p_w),color=w.sig.threshold),size=3)+
  geom_text(data=dataframeforplot1,aes(y=4.4+adjust,label=str_c('p.m=',p_m),color=m.sig.threshold),size=3)+
  geom_text(x=8,y=5.4+adjust,label=str_c('Num.Driver.Mutated = ',length(patient_cgc)),color="#E69F00",size=5)+
  geom_text(x=8,y=5.1+adjust,label=str_c('Num.Driver.Unmutated = ',length(patient_ncgc)),color="#999999",size=5)+
  scale_fill_manual(values=c( "#E69F00","#999999"))+
  scale_color_manual(values = c('#d11141','#f37735','#00aedb','#8c8c8c'),breaks=c('p<0.01','p<0.05','p<0.1','p>=0.1'))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 0,vjust = 1,size = 10),axis.title.y = element_text(size = 10))+
  labs(title='AUC Distributions Between Mutation Groups on NSCLC', y ="AUC", x = "Drug")
ggsave('cgc_AUC_box.pdf',width=12,height = 8)


mutation_p <- data.frame(value = drug_response_breast_all$lgec50, tag = drug_response_breast_all$tag, name = drug_response_breast_all$name)
mutaion_pvalues <- data.frame(drug=drug_names,p_t=1,p_w=1,p_m=1)
for (i in 1:length(drug_names)) {
  mutation_p1 <- mutation_p%>%subset(name==drug_names[i])
  mutation_p1$tag <- mutation_p1$tag%>%as.factor()
  t_test_value <- t.test(mutation_p1$value ~ mutation_p1$tag)
  mutaion_pvalues$p_t[i] <- t_test_value$p.value
  w_test_value <- wilcox.test(mutation_p1$value ~ mutation_p1$tag)
  mutaion_pvalues$p_w[i] <- w_test_value$p.value
  mutaion_pvalues$p_m[i] <- pvalue(median_test(mutation_p1$value ~ mutation_p1$tag))
}

dataframeforplot <- data.frame(value = drug_response_breast_all$lgec50, tag = drug_response_breast_all$tag, name = drug_response_breast_all$name)

# output_file <- str_c('box_',colnames(drug_response_breast_all)[i],'.pdf')


dataframeforplot1 <- dataframeforplot%>%distinct(name,tag)
dataframeforplot1$p_t <- round(mutaion_pvalues$p_t,3)
dataframeforplot1$p_w <- round(mutaion_pvalues$p_w,3)
dataframeforplot1$p_m <- round(mutaion_pvalues$p_m,3)
dataframeforplot1$sig.threshold <- as.factor(ifelse(dataframeforplot1$p_t<0.01,'p<0.01',ifelse(dataframeforplot1$p_t<0.05,'p<0.05',ifelse(dataframeforplot1$p_t<0.1,'p<0.1','p>=0.1'))))
dataframeforplot1$w.sig.threshold <- as.factor(ifelse(dataframeforplot1$p_w<0.01,'p<0.01',ifelse(dataframeforplot1$p_w<0.05,'p<0.05',ifelse(dataframeforplot1$p_w<0.1,'p<0.1','p>=0.1'))))
dataframeforplot1$m.sig.threshold <- as.factor(ifelse(dataframeforplot1$p_m<0.01,'p<0.01',ifelse(dataframeforplot1$p_m<0.05,'p<0.05',ifelse(dataframeforplot1$p_m<0.1,'p<0.1','p>=0.1'))))

dataframeforplot1$sig.threshold <- factor(dataframeforplot1$sig.threshold,levels =c('p<0.01','p<0.05','p<0.1','p>=0.1'))
dataframeforplot1$w.sig.threshold <- factor(dataframeforplot1$w.sig.threshold,levels =c('p<0.01','p<0.05','p<0.1','p>=0.1'))
dataframeforplot1$m.sig.threshold <- factor(dataframeforplot1$m.sig.threshold,levels =c('p<0.01','p<0.05','p<0.1','p>=0.1'))
adjust=0

dataframeforplot%>%ggplot(aes(x = name,y=value))+
  geom_boxplot(aes(fill = tag),color = 'black', alpha=0.4)+
  stat_summary(fun = mean,geom = "point", shape=5, size=1.5,aes(group=tag), position = position_dodge(0.8))+
  geom_text(data=dataframeforplot1,aes(y=4.8+adjust,label=str_c('p.t=',p_t),color=sig.threshold),size=3)+
  geom_text(data=dataframeforplot1,aes(y=4.6+adjust,label=str_c('p.w=',p_w),color=w.sig.threshold),size=3)+
  geom_text(data=dataframeforplot1,aes(y=4.4+adjust,label=str_c('p.m=',p_m),color=m.sig.threshold),size=3)+
  geom_text(x=8,y=5.4+adjust,label=str_c('Num.Driver.Mutated = ',length(patient_cgc)),color="#75FB53",size=5)+
  geom_text(x=8,y=5.1+adjust,label=str_c('Num.Driver.Unmutated = ',length(patient_ncgc)),color="#886fB9",size=5)+
  scale_fill_manual(values=c( "#75FB53", "#886fB9"))+
  scale_color_manual(values = c('#d11141','#f37735','#00aedb','#8c8c8c'),breaks=c('p<0.01','p<0.05','p<0.1','p>=0.1'))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 0,vjust = 1,size = 10),axis.title.y = element_text(size = 10))+
  labs(title='EC50 Distributions Between Mutation Groups on NSCLC', y ="log(EC50)", x = "Drug")
ggsave('cgc_EC50_box.pdf',width=12,height = 8)


if(!dir.exists('mutation_gene_drug')){
  dir.create('mutation_gene_drug')
}
setwd('mutation_gene_drug')


for (j in 1:length(cgc_breast$Gene.Symbol)) {
  patient_c <- rownames(mutation_matrix_breast)[mutation_matrix_breast[,cgc_breast$Gene.Symbol[j]]==1]
  patient_nc <- rownames(mutation_matrix_breast)[mutation_matrix_breast[,cgc_breast$Gene.Symbol[j]]==0]
  
  flag1 <- drug_response_breast$depmap_id%in%patient_c
  drug_response_breast_c <- drug_response_breast[flag1,]
  flag2 <- drug_response_breast$depmap_id%in%patient_nc
  drug_response_breast_nc <- drug_response_breast[flag2,]
  
  if(sum(flag1)>2&sum(flag2)>2){
    drug_response_breast_all1 <- drug_response_breast_c%>%rbind(drug_response_breast_nc)
    drug_response_breast_all1$tag <- c(rep(str_c(cgc_breast$Gene.Symbol[j],'_mutated'),nrow(drug_response_breast_c)),rep(str_c(cgc_breast$Gene.Symbol[j],'_unmutated'),nrow(drug_response_breast_nc)))
    
    # drug_response_breast_all1$tag <- c(rep('c',nrow(drug_response_breast_c)),rep('n',nrow(drug_response_breast_nc)))
    
    output_file1 <- str_c('box_auc',cgc_breast$Gene.Symbol[j],'.pdf')
    
    mutation_p <- data.frame(value = drug_response_breast_all1$auc_sd, tag = drug_response_breast_all1$tag, name = drug_response_breast_all1$name)
    mutaion_pvalues <- data.frame(drug=drug_names,p_t=1,p_w=1,p_m=1)
    for (i in 1:length(drug_names)) {
      mutation_p1 <- mutation_p%>%subset(name==drug_names[i])
      mutation_p1$tag <- mutation_p1$tag%>%as.factor()
      t_test_value <- t.test(mutation_p1$value ~ mutation_p1$tag)
      mutaion_pvalues$p_t[i] <- t_test_value$p.value
      w_test_value <- wilcox.test(mutation_p1$value ~ mutation_p1$tag)
      mutaion_pvalues$p_w[i] <- w_test_value$p.value
      mutaion_pvalues$p_m[i] <- pvalue(median_test(mutation_p1$value ~ mutation_p1$tag))
    }
    
    dataframeforplot <- data.frame(value = drug_response_breast_all1$auc_sd, tag = drug_response_breast_all1$tag, name = drug_response_breast_all1$name)
    
    dataframeforplot1 <- dataframeforplot%>%distinct(name,tag)
    dataframeforplot1$p_t <- round(mutaion_pvalues$p_t,3)
    dataframeforplot1$p_w <- round(mutaion_pvalues$p_w,3)
    dataframeforplot1$p_m <- round(mutaion_pvalues$p_m,3)
    dataframeforplot1$sig.threshold <- as.factor(ifelse(dataframeforplot1$p_t<0.01,'p<0.01',ifelse(dataframeforplot1$p_t<0.05,'p<0.05',ifelse(dataframeforplot1$p_t<0.1,'p<0.1','p>=0.1'))))
    dataframeforplot1$w.sig.threshold <- as.factor(ifelse(dataframeforplot1$p_w<0.01,'p<0.01',ifelse(dataframeforplot1$p_w<0.05,'p<0.05',ifelse(dataframeforplot1$p_w<0.1,'p<0.1','p>=0.1'))))
    dataframeforplot1$m.sig.threshold <- as.factor(ifelse(dataframeforplot1$p_m<0.01,'p<0.01',ifelse(dataframeforplot1$p_m<0.05,'p<0.05',ifelse(dataframeforplot1$p_m<0.1,'p<0.1','p>=0.1'))))
    
    dataframeforplot1$sig.threshold <- factor(dataframeforplot1$sig.threshold,levels =c('p<0.01','p<0.05','p<0.1','p>=0.1'))
    dataframeforplot1$w.sig.threshold <- factor(dataframeforplot1$w.sig.threshold,levels =c('p<0.01','p<0.05','p<0.1','p>=0.1'))
    dataframeforplot1$m.sig.threshold <- factor(dataframeforplot1$m.sig.threshold,levels =c('p<0.01','p<0.05','p<0.1','p>=0.1'))
    adjust=-1.5
    
    dataframeforplot%>%ggplot(aes(x = name,y=value))+
      geom_boxplot(aes(fill = tag),color = 'black', alpha=0.4)+
      stat_summary(fun = mean,geom = "point", shape=5, size=1.5,aes(group=tag), position = position_dodge(0.8))+
      geom_text(data=dataframeforplot1,aes(y=4.8+adjust,label=str_c('p.t=',p_t),color=sig.threshold),size=3)+
      geom_text(data=dataframeforplot1,aes(y=4.6+adjust,label=str_c('p.w=',p_w),color=w.sig.threshold),size=3)+
      geom_text(data=dataframeforplot1,aes(y=4.4+adjust,label=str_c('p.m=',p_m),color=m.sig.threshold),size=3)+
      geom_text(x=8,y=5.4+adjust,label=str_c('Num.',cgc_breast$Gene.Symbol[j],'.Mutated = ',length(patient_c)),color="#E69F00",size=5)+
      geom_text(x=8,y=5.1+adjust,label=str_c('Num.',cgc_breast$Gene.Symbol[j],'.Unmutated = ',length(patient_nc)),color="#999999",size=5)+
      scale_fill_manual(values=c( "#E69F00","#999999"))+
      scale_color_manual(values = c('#d11141','#f37735','#00aedb','#8c8c8c'),breaks=c('p<0.01','p<0.05','p<0.1','p>=0.1'))+
      theme_minimal()+
      theme(axis.text.x = element_text(angle = 0,vjust = 1,size = 10),axis.title.y = element_text(size = 10))+
      labs(title=str_c('AUC Distributions Between ',cgc_breast$Gene.Symbol[j], ' Mutation Groups on ',short[disease]),y ="AUC", x = "Drug")
    ggsave(output_file1,width=12,height = 8)
    
    # adjust=-1.5
    # 
    # dataframeforplot%>%ggplot(aes(x = name,y=value))+
    #   geom_boxplot(aes(fill = tag),color = 'black', alpha=0.4)+
    #   stat_summary(fun = mean,geom = "point", shape=5, size=1.5,aes(group=tag), position = position_dodge(0.8))+
    #   geom_text(data=dataframeforplot1,aes(y=4.8+adjust,label=str_c('p.t=',p_t),color=sig.threshold),size=3)+
    #   geom_text(data=dataframeforplot1,aes(y=4.6+adjust,label=str_c('p.w=',p_w),color=w.sig.threshold),size=3)+
    #   geom_text(data=dataframeforplot1,aes(y=4.4+adjust,label=str_c('p.m=',p_m),color=m.sig.threshold),size=3)+
    #   geom_text(x=8,y=5.4+adjust,label=str_c('Num.Driver.Mutated = ',length(patient_cgc)),color="#E69F00",size=5)+
    #   geom_text(x=8,y=5.1+adjust,label=str_c('Num.Driver.Unmutated = ',length(patient_ncgc)),color="#999999",size=5)+
    #   scale_fill_manual(values=c( "#E69F00","#999999"))+
    #   scale_color_manual(values = c('#d11141','#f37735','#00aedb','#8c8c8c'),breaks=c('p<0.01','p<0.05','p<0.1','p>=0.1'))+
    #   theme_minimal()+
    #   ylim(-3,4)+
    #   theme(axis.text.x = element_text(angle = 0,vjust = 1,size = 10),axis.title.y = element_text(size = 10))+
    #   labs(title='AUC Distributions Between Mutation Groups on NSCLC', y ="AUC", x = "Drug")
    # ggsave('cgc_AUC_box.pdf',width=12,height = 8)
    
    
    mutation_p <- data.frame(value = drug_response_breast_all1$lgec50, tag = drug_response_breast_all1$tag, name = drug_response_breast_all1$name)
    mutaion_pvalues <- data.frame(drug=drug_names,p_t=1,p_w=1,p_m=1)
    for (i in 1:length(drug_names)) {
      mutation_p1 <- mutation_p%>%subset(name==drug_names[i])
      mutation_p1$tag <- mutation_p1$tag%>%as.factor()
      t_test_value <- t.test(mutation_p1$value ~ mutation_p1$tag)
      mutaion_pvalues$p_t[i] <- t_test_value$p.value
      w_test_value <- wilcox.test(mutation_p1$value ~ mutation_p1$tag)
      mutaion_pvalues$p_w[i] <- w_test_value$p.value
      mutaion_pvalues$p_m[i] <- pvalue(median_test(mutation_p1$value ~ mutation_p1$tag))
    }
    
    dataframeforplot <- data.frame(value = drug_response_breast_all1$lgec50, tag = drug_response_breast_all1$tag, name = drug_response_breast_all1$name)
    
    dataframeforplot1 <- dataframeforplot%>%distinct(name,tag)
    dataframeforplot1$p_t <- round(mutaion_pvalues$p_t,3)
    dataframeforplot1$p_w <- round(mutaion_pvalues$p_w,3)
    dataframeforplot1$p_m <- round(mutaion_pvalues$p_m,3)
    dataframeforplot1$sig.threshold <- as.factor(ifelse(dataframeforplot1$p_t<0.01,'p<0.01',ifelse(dataframeforplot1$p_t<0.05,'p<0.05',ifelse(dataframeforplot1$p_t<0.1,'p<0.1','p>=0.1'))))
    dataframeforplot1$w.sig.threshold <- as.factor(ifelse(dataframeforplot1$p_w<0.01,'p<0.01',ifelse(dataframeforplot1$p_w<0.05,'p<0.05',ifelse(dataframeforplot1$p_w<0.1,'p<0.1','p>=0.1'))))
    dataframeforplot1$m.sig.threshold <- as.factor(ifelse(dataframeforplot1$p_m<0.01,'p<0.01',ifelse(dataframeforplot1$p_m<0.05,'p<0.05',ifelse(dataframeforplot1$p_m<0.1,'p<0.1','p>=0.1'))))
    
    
    dataframeforplot1$sig.threshold <- factor(dataframeforplot1$sig.threshold,levels =c('p<0.01','p<0.05','p<0.1','p>=0.1'))
    dataframeforplot1$w.sig.threshold <- factor(dataframeforplot1$w.sig.threshold,levels =c('p<0.01','p<0.05','p<0.1','p>=0.1'))
    dataframeforplot1$m.sig.threshold <- factor(dataframeforplot1$m.sig.threshold,levels =c('p<0.01','p<0.05','p<0.1','p>=0.1'))
    
    output_file2 <- str_c('box_ec50',cgc_breast$Gene.Symbol[j],'.pdf')
    adjust=0
    dataframeforplot%>%ggplot(aes(x = name,y=value))+
      geom_boxplot(aes(fill = tag),color = 'black', alpha=0.4)+
      stat_summary(fun = mean,geom = "point", shape=5, size=1.5,aes(group=tag), position = position_dodge(0.8))+
      geom_text(data=dataframeforplot1,aes(y=4.6+adjust,label=str_c('p.t=',p_t),color=sig.threshold),size=3)+
      geom_text(data=dataframeforplot1,aes(y=4.4+adjust,label=str_c('p.w=',p_w),color=w.sig.threshold),size=3)+
      geom_text(data=dataframeforplot1,aes(y=4.2+adjust,label=str_c('p.m=',p_m),color=m.sig.threshold),size=3)+
      geom_text(x=8,y=5.4+adjust,label=str_c('Num.',cgc_breast$Gene.Symbol[j],'.Mutated = ',length(patient_c)),color="#75FB53",size=5)+
      geom_text(x=8,y=5.1+adjust,label=str_c('Num.',cgc_breast$Gene.Symbol[j],'.Unmutated = ',length(patient_nc)),color="#886fB9",size=5)+
      scale_fill_manual(values=c( "#75FB53","#886fB9"))+
      scale_color_manual(values = c('#d11141','#f37735','#00aedb','#8c8c8c'),breaks=c('p<0.01','p<0.05','p<0.1','p>=0.1'))+
      theme_minimal()+
      theme(axis.text.x = element_text(angle = 0,vjust = 1,size = 10),axis.title.y = element_text(size = 10))+
      labs(title=str_c('EC50 Distributions Between ',cgc_breast$Gene.Symbol[j], ' Mutation Groups on ',short[disease]), y ="log(EC50)", x = "Drug")
    
    ggsave(output_file2,width=12,height = 8)
  }
}

setwd('..')


# cnv
if(!dir.exists('cnv')){
  dir.create('cnv')
}
setwd('cnv')

flag <- c()
for (gene in cgc_breast$Gene.Symbol) {
  temp <- colnames(ccle_cnv)%>%str_which(str_c(gene,'\\.'))
  flag <- flag%>%append(temp[1])
}
cgc_breast <- cgc_breast[!is.na(flag),]

auc_cor_table <- data_frame(cancer=NULL,drug=NULL,gene=NULL, pearson=NULL, spearman=NULL,kendall=NULL)
ec50_cor_table <- data_frame(cancer=NULL,drug=NULL,gene=NULL, pearson=NULL, spearman=NULL,kendall=NULL)
auc_fit_table <- data_frame(cancer=NULL,drug=NULL,gene=NULL, lm_result_aic=NULL,gam_result_aic=NULL,rlm_result_aic=NULL,
                            lm_result_r2=NULL,gam_result_r2=NULL,rlm_result_mse=NULL)
ec50_fit_table <- data_frame(cancer=NULL,drug=NULL,gene=NULL, lm_result_aic=NULL,gam_result_aic=NULL,rlm_result_aic=NULL,
                             lm_result_r2=NULL,gam_result_r2=NULL,rlm_result_mse=NULL)

for (j in 1:length(drug_names)) {
  if(!dir.exists(drug_names[j])){
    dir.create(drug_names[j])
  }
  setwd(drug_names[j])
  
  drug_response_breast_drug <- drug_response_breast%>%subset(name==drug_names[j])
  id_list <- ccle_cnv$X%>%intersect(drug_response_breast_drug$depmap_id)
  drug_response_breast_drug <- drug_response_breast_drug[drug_response_breast_drug$depmap_id%in%id_list,]
  ccle_cnv_breast <- ccle_cnv[ccle_cnv$X%in%id_list,]
  
  flag <- drug_response_breast_drug$depmap_id%>%match(ccle_cnv_breast$X)
  ccle_cnv_breast <- ccle_cnv_breast[flag,]
  
  flag <- c()
  for (gene in cgc_breast$Gene.Symbol) {
    temp <- colnames(ccle_cnv)%>%str_which(str_c(gene,'\\.'))
    flag <- flag%>%append(temp[1])
  }
  flag <- flag[!is.na(flag)]
  ccle_cnv_breast <- ccle_cnv_breast[,flag]
  
  ccle_cnv_breast$ec50 <- drug_response_breast_drug$lgec50
  ccle_cnv_breast$auc <- drug_response_breast_drug$auc_sd
  
  
  for (i in 1:length(cgc_breast$Gene.Symbol)) {
    
    cor_pearson <- cor.test(as.vector(ccle_cnv_breast[,i]),as.vector(ccle_cnv_breast$auc),method = 'pearson')
    cor_spearman <- cor.test(as.vector(ccle_cnv_breast[,i]),as.vector(ccle_cnv_breast$auc),method = 'spearman')
    cor_kendall <- cor.test(as.vector(ccle_cnv_breast[,i]),as.vector(ccle_cnv_breast$auc),method = 'kendall')
    auc_cor_table1 <- data_frame(cancer=short[disease],drug=drug_names[j],gene=cgc_breast$Gene.Symbol[i], pearson=cor_pearson$estimate%>%round(3),
                                 spearman=cor_spearman$estimate%>%round(3),kendall=cor_kendall$estimate%>%round(3))
    auc_cor_table <- auc_cor_table%>%rbind(auc_cor_table1)
    auc_cor_table <- auc_cor_table%>%arrange(desc(abs(pearson)))
    
    lm_result <- lm(as.vector(ccle_cnv_breast[,i]) ~ as.vector(ccle_cnv_breast$auc))
    rlm_result <- rlm(as.vector(ccle_cnv_breast[,i]) ~ as.vector(ccle_cnv_breast$auc))
    gam_result <- gam(as.vector(ccle_cnv_breast[,i]) ~ s(as.vector(ccle_cnv_breast$auc),bs='cr'))
    
    auc_fit_table1 <- data_frame(cancer=short[disease],drug=drug_names[j],gene=cgc_breast$Gene.Symbol[i], lm_result_aic=AIC(lm_result),gam_result_aic=AIC(gam_result),rlm_result_aic=AIC(rlm_result),
                                 lm_result_r2=summary(lm_result)$r.sq,gam_result_r2=summary(gam_result)$r.sq,rlm_result_mse=rlm_result$s)
    auc_fit_table <- auc_fit_table%>%rbind(auc_fit_table1)
    
    tab <- as.data.frame(
      c(
        P = cor_pearson$estimate%>%round(3),
        S = cor_spearman$estimate%>%round(3),
        K = cor_kendall$estimate%>%round(3)
      )
    )
    
    output_file <- str_c('aucpoint_',drug_names[j],colnames(ccle_cnv_breast)[i],'.pdf')
    dataframeforplot <- data.frame(cnv = ccle_cnv_breast[,i], auc = ccle_cnv_breast$auc)
    p <- ggplot(dataframeforplot, aes(x=cnv,y=auc))+
      geom_point(color="#00AFBB", alpha=0.6,size=3)+
      geom_smooth(method=lm,se = T,aes(color="red"),fill="red",alpha=0.13)+
      geom_smooth(method='rlm',formula=y ~ splines::bs(x,3),se = T,aes(color="turquoise"),fill="turquoise",alpha=0.13)+
      geom_smooth(method='gam',formula = y ~ rms::rcs(x, 3),se = T,aes(color="orange"),fill="orange",alpha=0.13)+
      theme_minimal()+
      scale_color_manual(name='Fitting',values=c("orange","red","turquoise"), labels=c('GAM','LM','RLM'),)+
      # scale_fill_manual(name='Fitting',values=c("red","turquoise","orange"), labels=c('lm','glm','gam'),)+
      labs(title=str_c(cgc_breast$Gene.Symbol[i],' CNV & ',drug_names[j],' AUC '),y ="AUC", x = "CNV")+
      theme(legend.position = c(1.1,0.35))
    
    
    ggdraw()+
      draw_plot(p,width = 0.8)+
      draw_plot(tableGrob(unname(tab)),x=0.8,y=0.14, width = 0.15)
    
    ggsave(output_file,width=6,height = 5)
    
    
    cor_pearson <- cor.test(as.vector(ccle_cnv_breast[,i]),as.vector(ccle_cnv_breast$ec50),method = 'pearson')
    cor_spearman <- cor.test(as.vector(ccle_cnv_breast[,i]),as.vector(ccle_cnv_breast$ec50),method = 'spearman')
    cor_kendall <- cor.test(as.vector(ccle_cnv_breast[,i]),as.vector(ccle_cnv_breast$ec50),method = 'kendall')
    ec50_cor_table1 <- data_frame(cancer=short[disease],drug=drug_names[j],gene=cgc_breast$Gene.Symbol[i], pearson=cor_pearson$estimate%>%round(3),
                                  spearman=cor_spearman$estimate%>%round(3),kendall=cor_kendall$estimate%>%round(3))
    ec50_cor_table <- ec50_cor_table%>%rbind(ec50_cor_table1)
    ec50_cor_table <- ec50_cor_table%>%arrange(desc(abs(pearson)))
    
    lm_result <- lm(as.vector(ccle_cnv_breast[,i]) ~ as.vector(ccle_cnv_breast$ec50))
    rlm_result <- rlm(as.vector(ccle_cnv_breast[,i]) ~ as.vector(ccle_cnv_breast$ec50))
    gam_result <- gam(as.vector(ccle_cnv_breast[,i]) ~ s(as.vector(ccle_cnv_breast$ec50),bs='cr'))
    
    ec50_fit_table1 <- data_frame(cancer=short[disease],drug=drug_names[j],gene=cgc_breast$Gene.Symbol[i], lm_result_aic=AIC(lm_result),gam_result_aic=AIC(gam_result),rlm_result_aic=AIC(rlm_result),
                                  lm_result_r2=summary(lm_result)$r.sq,gam_result_r2=summary(gam_result)$r.sq,rlm_result_mse=rlm_result$s)
    ec50_fit_table <- ec50_fit_table%>%rbind(ec50_fit_table1)
    
    tab <- as.data.frame(
      c(
        P = cor_pearson$estimate%>%round(3),
        S = cor_spearman$estimate%>%round(3),
        K = cor_kendall$estimate%>%round(3)
      )
    )
    
    output_file <- str_c('ec50point_',drug_names[j],colnames(ccle_cnv_breast)[i],'.pdf')
    dataframeforplot <- data.frame(cnv = ccle_cnv_breast[,i], EC50 = ccle_cnv_breast$ec50)
    temp <- dataframeforplot$EC50%>%quantile()
    sup <-temp[4]+1.5*(temp[4]-temp[2])
    dataframeforplot <- dataframeforplot%>%subset(EC50<sup)
    p <- ggplot(dataframeforplot, aes(x=cnv,y=EC50))+
      geom_point(color="#E7B800", alpha=0.6,size=3)+
      geom_smooth(method=lm,se = T,aes(color="red"),fill="red",alpha=0.13)+
      geom_smooth(method='rlm',formula=y ~ splines::bs(x,3),se = T,aes(color="turquoise"),fill="turquoise",alpha=0.13)+
      geom_smooth(method='gam',formula = y ~ rms::rcs(x, 3),se = T,aes(color="orange"),fill="orange",alpha=0.13)+
      theme_minimal()+
      scale_color_manual(name='Fitting',values=c("orange","red","turquoise"), labels=c('GAM','LM','RLM'),)+
      # scale_fill_manual(name='Fitting',values=c("red","turquoise","orange"), labels=c('lm','glm','gam'),)+
      labs(title=str_c(cgc_breast$Gene.Symbol[i],' CNV & ',drug_names[j],' EC50 '),y ="log(EC50)", x = "CNV")+
      theme(legend.position = c(1.1,0.35))
    
    ggdraw()+
      draw_plot(p,width = 0.8)+
      draw_plot(tableGrob(unname(tab)),x=0.8,y=0.14, width = 0.15)
    
    ggsave(output_file,width=6,height = 5)
  }
  
  setwd('..')
}
write.csv(auc_cor_table,'auc_cor_table.csv')
write.csv(ec50_cor_table,'ec50_cor_table.csv')
write.csv(auc_fit_table,'auc_fit_table.csv')
write.csv(ec50_fit_table,'ec50_fit_table.csv')



# if(!dir.exists('cnv_density')){
#   dir.create('cnv_density')
# }
# setwd('cnv_density')

# drug_response_breast <- drug_response%>%subset(indication=='breast cancer')
id_list <- ccle_cnv$X%>%intersect(drug_response_breast$depmap_id)
drug_response_breast_sub <- drug_response_breast[drug_response_breast$depmap_id%in%id_list,]
ccle_cnv_breast <- ccle_cnv[ccle_cnv$X%in%id_list,]

flag <- drug_response_breast_sub$depmap_id%>%match(ccle_cnv_breast$X)
ccle_cnv_breast <- ccle_cnv_breast[flag,]

flag <- c()
for (gene in cgc_breast$Gene.Symbol) {
  temp <- colnames(ccle_cnv)%>%str_which(str_c(gene,'\\.'))
  flag <- flag%>%append(temp[1])
}
flag <- flag[!is.na(flag)]
ccle_cnv_breast <- ccle_cnv_breast[,flag]

ccle_cnv_breast$lgec50 <- drug_response_breast_sub$lgec50
ccle_cnv_breast$auc <- drug_response_breast_sub$auc_sd
ccle_cnv_breast$drug <- drug_response_breast_sub$name

# for (i in 1:length(cgc_breast$Gene.Symbol)) {
#   dataframeforplot <- data.frame(gene = ccle_cnv_breast[,i])
#   output_file <- str_c(cgc_breast$Gene.Symbol[i],'.pdf')
#   dataframeforplot%>%ggplot(aes(x=gene,..density../length(drug_names)))+
#     geom_histogram(colour="black", fill="white")+
#     geom_density(alpha=0.2, fill = "#E69F00")+
#     theme_minimal()+
#     labs(title=cgc_breast$Gene.Symbol[i],x='CNV',y='Density')+
#     theme(plot.title = element_text(hjust=0.5))
#   
#   ggsave(output_file,width=3,height = 3)
# }
# setwd('..')

if(!dir.exists('cnv_box')){
  dir.create('cnv_box')
}
setwd('cnv_box')

ccle_cnv_breast$cnv_max <- 0
for (i in 1:nrow(ccle_cnv_breast)) {
  ccle_cnv_breast$cnv_max[i] <- round(max(ccle_cnv_breast[i,1:length(cgc_breast$Gene.Symbol)]))
}
ccle_cnv_breast$cnv_max[ccle_cnv_breast$cnv_max>=3] <- 3
ccle_cnv_breast$cnv_max <- ccle_cnv_breast$cnv_max%>%as.character()
ccle_cnv_breast$cnv_max[ccle_cnv_breast$cnv_max=='3'] <- '3+'
drug_names <- unique(ccle_cnv_breast$drug)
for (j in 1:length(drug_names)) {
  ccle_cnv_breast_sub <- ccle_cnv_breast%>%subset(drug==drug_names[j])
  
  dataframeforplot1 <- ccle_cnv_breast_sub$cnv_max%>%table()%>%as.data.frame()
  colnames(dataframeforplot1) <- c('cnv_max','num')
  
  output_file <- str_c('auc_',drug_names[j],'.pdf')
  ccle_cnv_breast_sub%>%ggplot(aes(x=cnv_max,y=auc))+
    geom_boxplot(aes(fill = cnv_max))+
    geom_text(data=dataframeforplot1,aes(y=Inf,label=str_c('n = ',num),vjust=2),size=4,color='#f37735')+
    theme_minimal()+
    scale_fill_manual(values=wes_palette(n=4, name="GrandBudapest2"))+
    labs(title=drug_names[j], y ="AUC", x = "CNV")+
    theme(plot.title = element_text(hjust=0.5))
  ggsave(output_file,width=4,height = 3)
  
  auc_med <- ccle_cnv_breast_sub%>%group_by(cnv_max)%>%summarise(value=median(auc),ymin=quantile(auc,probs = 0.25),ymax=quantile(auc,probs = 0.75))
  
  f <- ggplot(auc_med, aes(x = cnv_max, y = value,
                           ymin = ymin, ymax = ymax))
  
  f + geom_line(aes(group=1)) + 
    geom_point()+ 
    geom_errorbar(width=.2,position = position_dodge(0.05))+
    theme_minimal()+
    labs(title=drug_names[j], y ="AUC", x = "CNV")+
    theme(plot.title = element_text(hjust=0.5))
  ggsave(str_c('line_auc_',drug_names[j],'.pdf'),width=4,height = 3)
  
  output_file <- str_c('ec50_',drug_names[j],'.pdf')
  ccle_cnv_breast_sub%>%ggplot(aes(x=cnv_max,y=lgec50))+
    geom_boxplot(aes(fill = cnv_max))+
    geom_text(data=dataframeforplot1,aes(y=Inf,label=str_c('n = ',num),vjust=2),size=4,color='#f37735')+
    theme_minimal()+
    scale_fill_manual(values=wes_palette(n=4, name="GrandBudapest2"))+
    labs(title=drug_names[j], y ="log(EC50)", x = "CNV")+
    theme(plot.title = element_text(hjust=0.5))
  ggsave(output_file,width=4,height = 3)
  
  auc_med <- ccle_cnv_breast_sub%>%group_by(cnv_max)%>%summarise(value=median(lgec50),ymin=quantile(lgec50,probs = 0.25),ymax=quantile(lgec50,probs = 0.75))
  
  f <- ggplot(auc_med, aes(x = cnv_max, y = value,
                           ymin = ymin, ymax = ymax))
  
  f + geom_line(aes(group=1)) + 
    geom_point()+ 
    geom_errorbar(width=.2,position = position_dodge(0.05))+
    theme_minimal()+
    labs(title=drug_names[j], y ="log(EC50)", x = "CNV")+
    theme(plot.title = element_text(hjust=0.5))
  ggsave(str_c('line_ec50_',drug_names[j],'.pdf'),width=4,height = 3)
  
}

setwd('..')
