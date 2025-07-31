rm(list = ls());gc()
library(survival)
library(randomForestSRC)
library(glmnet)
library(plsRcox)
library(superpc)
library(gbm)
library(CoxBoost)
library(survivalsvm)
library(dplyr)
library(tibble)
library(BART)
library(limma)
library(tidyverse)
library(dplyr)
library(miscTools)
library(compareC)
library(ggplot2)
library(ggsci)
library(tidyr)
library(ggbreak)
library(data.table)
library(Mime1)
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 
setwd('D:/SynologyDrive/R/LGG/')
#---数据读取+取交集----
sets <- c('TCGA',
          'CGGA325',
          'CGGA301',
          'CGGA693',
          'GSE16011',
          'MTAB3892')
datasets <- list()
input_sets <- lapply(1:6, function(x){
  exp <- readRDS(file = paste0('cleandata/',sets[x],'_exp.rds'))
  clin <- readRDS(file = paste0('cleandata/',sets[x],'_clin.rds')) 
  df <- data.frame(
    ID=rownames(clin),
    OS.time=clin$OS.time,
    OS=clin$OS
  )
  
  df <- cbind(df,t(exp))
  df <- df[df$OS.time>0,]
  datasets[[x]] <- as.data.frame(df)
}
)
names(input_sets) <- c('TCGA',
                       'CGGA325',
                       'CGGA301',
                       'CGGA693',
                       'GSE16011',
                       'MTAB3892')
genes <- clusterProfiler::read.gmt('data/WP_ELECTRON_TRANSPORT_CHAIN_OXPHOS_SYSTEM_IN_MITOCHONDRIA.v2024.1.Hs.gmt')
genelist <- tinyarray::intersect_all(
  colnames(input_sets[[1]]),colnames(input_sets[[2]]),colnames(input_sets[[3]]),
  colnames(input_sets[[4]]),colnames(input_sets[[5]]),colnames(input_sets[[6]]),
  genes$gene_symbol
)
new_input_sets <- lapply(1:6, function(x){
  input_sets[[x]][,c(1:2)] <- apply(input_sets[[x]][,c(1:2)],2,as.factor)
  input_sets[[x]][,c(2:ncol(input_sets[[x]]))] <- apply(input_sets[[x]][,c(2:ncol(input_sets[[x]]))],2,as.numeric)
  input_sets[[x]] <- as.data.frame(input_sets[[x]])
})
names(new_input_sets) <- c('TCGA',
                           'CGGA325',
                           'CGGA301',
                           'CGGA693',
                           'GSE16011',
                           'MTAB3892')
saveRDS(new_input_sets,file = 'model_input.rds')
res <- ML.Dev.Prog.Sig(train_data = new_input_sets$TCGA,
                       list_train_vali_Data = new_input_sets,
                       unicox.filter.for.candi = T,
                       unicox_p_cutoff = 0.001,
                       candidate_genes = genelist,
                       mode = 'all',nodesize =5,
                       seed = 5201314)
saveRDS(res,file = 'res_ML.rds')
colnames(res$ml.res$`StepCox[forward] + plsRcox`$dataX)
source('cindex_dis_all2.R')
cindex_dis_all2(res,validate_set = names(new_input_sets)[-1],
               dataset_col = c("#E64B35C6", "#ED00007F", "#CD202C7F", "#C500847F", "#7D5CC67F", "#0073C27F"),
               order =names(new_input_sets),width = 0.15)
cowplot::ggsave2(filename = 'plot/Model_cindex.pdf',height = 20,width = 8)
source('cindex_dis_select2.R')
cindex_dis_select2(res,dataset_col = c("#E64B35C6", "#ED00007F", "#CD202C7F", "#C500847F", "#7D5CC67F", "#0073C27F"),
                  model="StepCox[forward] + plsRcox",
                  order= names(new_input_sets))
cowplot::ggsave2(filename = 'plot/Model_cindex_Fig1.pdf',height = 6,width = 8)

vars <- colnames(res$ml.res$`StepCox[forward] + plsRcox`$dataX)
ggsci::pal_npg('nrc',alpha = 0.78)(2)
survplot <- vector("list",6) 
for (i in c(1:6)) {
  print(survplot[[i]]<-rs_sur(res, model_name = "StepCox[forward] + plsRcox",
                              dataset = names(new_input_sets)[i],
                              color=c("#E64B35C6","#4DBBD5C6"),
                              median.line = "hv",
                              cutoff = 0.5,
                              conf.int = T,
                              xlab="Days",pval.coord=c(100,0.2)))
}
pdf(file = 'plot/surv_model.pdf',height = 8,width = 12)
aplot::plot_list(gglist=survplot,ncol=3)
dev.off()
#### 特定队列的time-dependent的auc ####

lapply(1:6, function(x){
  all.auc.1y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = new_input_sets[[x]],
                               inputmatrix.list = new_input_sets,mode = 'all',AUC_time = 1,
                               auc_cal_method="KM")
  all.auc.3y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = new_input_sets[[x]],
                               inputmatrix.list = new_input_sets,mode = 'all',AUC_time = 3,
                               auc_cal_method="KM")
  all.auc.5y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = new_input_sets[[x]],
                               inputmatrix.list = new_input_sets,mode = 'all',AUC_time = 5,
                               auc_cal_method="KM")
  auc_dis_select(list(all.auc.1y,all.auc.3y,all.auc.5y),
                 model_name="StepCox[forward] + plsRcox",
                 dataset = names(new_input_sets),
                 order= names(new_input_sets),
                 year=c(1,3,5))
  cowplot::ggsave2(filename = paste0('plot/tdAUC_',names(new_input_sets[x]),'.pdf'),height = 5,width = 5)
})
### AUC图###
library(timeROC)
library(survival)
for (i in 1:6) {
  new_dat <- as.data.frame(res$riskscore$`StepCox[forward] + plsRcox`[[i]])
  result <-with(new_dat, timeROC(T=OS.time,
                                  delta=OS,
                                  marker=RS,
                                  cause=1,
                                  times=c(365,1095,1825),
                                  iid = TRUE))
  dat = data.frame(fpr = as.numeric(result$FP),
                   tpr = as.numeric(result$TP),
                   time = rep(as.factor(c(365,1095,1825)),each = nrow(result$TP)))
  ggplot() +
    geom_line(data = dat,aes(x = fpr, y = tpr,color = time),size = 1) +
    scale_color_manual(name = NULL,values = c("#7AA6DCCC", "#A73030CC", "#003C67CC"),
                       labels = paste0("AUC of ",c(365/365,1095/365,1825/365),"-y survival: ",
                                       format(round(result$AUC,2),nsmall = 2)))+
    geom_line(aes(x=c(0,1),y=c(0,1)),color = "grey")+
    theme_bw()+
    theme(panel.grid = element_blank(),
          axis.title = element_text(size=20),
          axis.text = element_text(size=18,colour = 'black'),
          legend.text = element_text(size = 24),
          legend.background = element_rect(linetype = 1, size = 0.2, colour = "black"),
          legend.position = c(0.695,0.125))+
    scale_x_continuous(expand = c(0.005,0.005))+
    scale_y_continuous(expand = c(0.005,0.005))+
    labs(x = "1 - Specificity",
         y = "Sensitivity")+
    coord_fixed()
  ggsave(last_plot(),filename = paste0('plot/tdAUC_',names(res$riskscore$`StepCox[forward] + plsRcox`)[i],'.pdf'),
         width = 8,height = 8)
  
}


####特定模型单因素cox回归的Meta分析####
unicox.rs.res <- cal_unicox_ml_res(res.by.ML.Dev.Prog.Sig = res,
                                   optimal.model = "StepCox[forward] + plsRcox",type ='categorical')
metamodel <- cal_unicox_meta_ml_res(input = unicox.rs.res)
pdf('plot/Meta_unicox_fig1.pdf',height = 5,width = 10)
meta_unicox_vis(metamodel,dataset_col = c("#E64B35C6", "#ED00007F", "#CD202C7F", "#C500847F", "#7D5CC67F", "#0073C27F"),
                dataset = names(new_input_sets))
dev.off()
###将特定模型的HR与先前发布的模型进行比较####
source('HR_com2.R')
rs.glioma.lgg.gbm <- cal_RS_pre.prog.sig(use_your_own_collected_sig = F,
                                         type.sig = c('LGG','Glioma'),
                                         list_input_data = new_input_sets)
HR_com2(rs.glioma.lgg.gbm,color = c("#0073C2CC",'white',"#A73030CC"),
       res,dataset_col = c("#E64B35C6", "#ED00007F", "#CD202C7F", "#C500847F", "#7D5CC67F", "#0073C27F"),
       model_name="StepCox[forward] + plsRcox",
       dataset=names(new_input_sets),
       type = "categorical")
cowplot::ggsave2(filename = 'plot/Model_index_compare_Fig2.pdf',height = 2,width = 10)

###将特定模型的C 指数与先前发布的模型进行比较  Figure 2B####

cc.glioma.lgg.gbm <- cal_cindex_pre.prog.sig(use_your_own_collected_sig = F,
                                             type.sig = c('LGG','Glioma'),
                                             list_input_data = new_input_sets)
cindex_comp(cc.glioma.lgg.gbm,
            res,dataset_col = c("#E64B35C6", "#ED00007F", "#CD202C7F", "#C500847F", "#7D5CC67F", "#0073C27F"),
            model_name="StepCox[forward] + plsRcox",
            dataset=names(new_input_sets))
cowplot::ggsave2(filename = 'plot/Model_Cindex_compare_Fig2.pdf',height = 10,width = 15)
str(new_input_sets[['Dataset6']][1:3])

#### Clinical C-index Figure 2C #####
library(data.table)
library(ggh4x)
source('clin_cox.R')
new_input_sets <- readRDS('model_input.rds')
res <- readRDS('res_ML.rds')
sets <- c('TCGA','CGGA325','CGGA301','CGGA693','GSE16011','MTAB3892') #输入需要纳入得数据集
clin_vars <- c('Age','Gender','Grade','IDH','PQ') #输入需要纳入的变量名

clin_sets <- list()
clin_sets <- lapply(1:length(sets), function(x){
  clin <- readRDS(file = paste0('cleandata/',sets[x],'_clin.rds'))
  clin <- clin %>% 
      dplyr::mutate(Age=case_when(Age > 40 ~ 1,
                              Age <= 40 ~ 0))
  df2 <- new_input_sets[[x]][1:3]
  clin <- clin[match(df2$ID,rownames(clin)),]
  df <- cbind(df2,RS=res$riskscore$`StepCox[forward] + plsRcox`[[x]][,4],clin[,clin_vars])
  clin_sets[[x]] <- df
})
names(clin_sets) <- sets
vars <- c('RS','Age','Gender','Grade','IDH','PQ') #输入需要纳入的变量名
result <- cal_index_clin(data = clin_sets,
                       sets = sets,
                       vars = vars)
dataset_col <- c("#E64B35C6", "#ED00007F", "#CD202C7F", "#C500847F", "#7D5CC67F", "#0073C27F")
result$Cindex <- round(result$Cindex,2)
result$Datasets <- factor(result$Datasets,levels = sets)

p1 <- result %>% 
  ggplot( aes(x = Variables, y = Cindex, fill = Datasets)) +
  geom_bar(position = "dodge", stat = "identity",color = "black") +
  
  geom_errorbar(aes(ymin = Cindex - se, ymax = Cindex + se),
                position = position_dodge(0.9), width = .2)+
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.3),
    legend.position = "",
    plot.title = element_text(hjust = 0.5),
    axis.text.y = element_text(size = 15,colour = 'black',face = 'plain'),
    panel.background = element_rect(fill = "white")
  ) +
  labs(y = "Cindex", x = "") +
  coord_flip()
ridiculous_strips <- strip_themed(text_x = elem_list_text(size=12),
  background_x = elem_list_rect(
    fill = 'white'),
  by_layer_x = F)
rs_index <- data.frame(Datasets=factor(levels(result$Datasets)),
                       Cindex=result$Cindex[result$Variables=='RS']) 
  
p1 +   
  facet_grid2( ~ Datasets, strip = ridiculous_strips)+
  scale_fill_manual(values = dataset_col, name = 'Datasets')+
  geom_hline(aes(yintercept = Cindex), rs_index,linetype="dashed",
             size=0.5, colour="gray56")
cowplot::ggsave2(filename = 'plot/Clin_cindex_datasets.pdf',height = 6,width = 12)


### cibersort ####
library(ggplot2)
library(tidyverse)
library(ggsci)
library(ggpubr)
res <- readRDS('res_ML.rds')
devo <- readRDS('devo.rds')
gene <- read.table('immunomodulator.txt',sep = '\t')
sets <- c('TCGA','CGGA325','CGGA301','CGGA693','GSE16011','MTAB3892')
for(j in sets){
  dat <- devo[[j]]$tme_combine$cibersort
  RS <- res$riskscore$`StepCox[forward] + plsRcox`[[j]]$RS
  names(RS) <- res$riskscore$`StepCox[forward] + plsRcox`[[j]]$ID
  group <- factor(ifelse(RS>median(RS),'High','Low'),levels=c('High','Low'))
  data2 <- gather(dat,key = 'samples',value = Score,colnames(dat)[2:ncol(dat)]) 
  group <- group[match(data2$samples,names(group),nomatch = F)]
  data2$group <- group
  
  ggplot(as.data.frame(data2),aes(cell_type,Score,fill = group)) + 
    geom_boxplot(outlier.shape = NA) + 
    theme_bw() +
    labs(x = "", y = "Fraction") +
    theme(legend.position = "top") + 
    theme(axis.text.x = element_text(size=14,colour = 'black'),
          axis.text.y = element_text(size=14,colour = 'black'),
          legend.text = element_text(size=14),
          legend.title = element_text(size=16),
          axis.title.y = element_text(size=16),
          axis.title.x = element_text(size=20),
          axis.line = element_line(size = 0.5))+
    scale_fill_npg(alpha = 0.75)+ 
    stat_compare_means(aes(group = group,label = ..p.signif..),method = "wilcox.test")+
    coord_flip()
  cowplot::ggsave2(filename = paste0('plot/LGG_Cib_',j,'.pdf'),width = 6,height = 12)
  dev.off()
  
}

result <- do.call(cbind, lapply(names(devo), function(dataset) {
  # 检查 tme_combine 和 cibersort 是否存在
  if (!is.null(devo[[dataset]]$tme_combine) &&
      "cibersort" %in% names(devo[[dataset]]$tme_combine)) {
    cibersort_data <- devo[[dataset]]$tme_combine$cibersort
    
    # 如果不是第一个数据集，剔除重复的 "cell type" 列
    if ("cell_type" %in% colnames(cibersort_data)) {
      cibersort_data <- cibersort_data[, !(colnames(cibersort_data) == "cell_type"), drop = FALSE]
    }
    return(cibersort_data)
  } else {
    return(NULL) # 忽略不存在的数据集
  }
}))

result2 <- do.call(rbind, lapply(names(res$riskscore$`StepCox[forward] + plsRcox`), function(dataset) {
  if (!is.null(res$riskscore$`StepCox[forward] + plsRcox`[[dataset]]) &&
      "RS" %in% names(res$riskscore$`StepCox[forward] + plsRcox`[[dataset]])) {
    return(
      data.frame(RS=res$riskscore$`StepCox[forward] + plsRcox`[[dataset]]$RS,
                 ID=res$riskscore$`StepCox[forward] + plsRcox`[[dataset]]$ID,
                 dataset=dataset)
           )
  } else {
    return(NULL) # 忽略不存在的数据集
  }
}))
result2 <- result2 %>% 
  filter(!duplicated(ID)) %>% 
  group_by(dataset) %>% 
  mutate(group=ifelse(RS>median(RS),'High','Low')) 
 
# 查看合并后的数据框
print(result)

# 查看合并后的数据框
colnames(result)

result <- cbind(cell_type=devo$TCGA$tme_combine$cibersort$cell_type,result)
result <- result[!duplicated(colnames(result)),]
result <- result[1:22,]
data2 <- gather(result,key = 'samples',value = Score,colnames(result)[2:ncol(result)]) %>% 
  select(cell_type,samples,Score)

groups <- result2[match(data2$samples,result2$ID,nomatch = F),]
data2$group <- groups$group


ggplot(as.data.frame(data2),aes(cell_type,Score,fill = group)) + 
  geom_boxplot(outlier.shape = NA) + 
  theme_bw() +
  labs(x = "", y = "Fraction") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(size=14,colour = 'black'),
        axis.text.y = element_text(size=14,colour = 'black'),
        legend.text = element_text(size=14),
        legend.title = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.title.x = element_text(size=20),
        axis.line = element_line(size = 0.5))+
  scale_fill_npg(alpha = 0.75)+ 
  stat_compare_means(aes(group = group,label = ..p.signif..),method = "wilcox.test")+
  coord_flip()
cowplot::ggsave2(filename = paste0('plot/LGG_Cib_mix.pdf'),width = 6,height = 12)
dev.off()


### ic50 #####
library(oncoPredict)
library(parallel)
library(ggpubr)
library(reshape2)
set.seed(999)

GDSC2_Expr =readRDS(file='./DataFiles/Training Data/GDSC2_Expr (RMA Normalized and Log Transformed).rds')
GDSC2_Res=readRDS(file='./DataFiles/Training Data/GDSC2_Res.rds')
GDSC2_Res <- exp(GDSC2_Res)
ic50_GDSC2 <- list()
ic50_GDSC2 <- list()
for(j in sets){
  dat <- new_input_sets[[j]]
  RS <- res$riskscore$`StepCox[forward] + plsRcox`[[j]]$RS
  names(RS) <- res$riskscore$`StepCox[forward] + plsRcox`[[j]]$ID
  group <- factor(ifelse(RS>median(RS),'High','Low'),levels=c('Low','High'))
  group <- group[match(dat$ID,names(group))]
  rownames(dat) <- NULL
  rownames(dat) <- names(group)
  dat <- as.data.frame(t(dat[,4:ncol(dat)]))
  testExpr <- as.matrix(dat)
  if(! dir.exists(paste0('D:/SynologyDrive/R/LY/LGG/',j))){
    dir.create(paste0('D:/SynologyDrive/R/LGG/',j))
  }
  setwd(paste0('D:/SynologyDrive/R/LGG/',j))
  calcPhenotype(trainingExprData = GDSC2_Expr,
                trainingPtype = GDSC2_Res,
                testExprData = testExpr,
                batchCorrect = 'eb',  #                 
                powerTransformPhenotype = TRUE,              
                removeLowVaryingGenes = 0.2,              
                minNumSamples = 10,               
                printOutput = TRUE,               
                removeLowVaringGenesFrom = 'rawData')
  ic50_GDSC2[[j]] <- read.csv('./calcPhenotype_Output/DrugPredictions.csv', header = T ,stringsAsFactors = F ,check.names = F)
  names(ic50_GDSC2[j]) <- j
}
saveRDS(ic50_GDSC2,file = 'ic50.rds')
#Drug selection
  setwd('D:/SynologyDrive/R/LGG/')
  drug_names <- colnames(read.csv(paste0('MTAB3892/calcPhenotype_Output/DrugPredictions.csv'), header = T ,stringsAsFactors = F ,check.names = F))
  drug_names <- drug_names[-1]
  ic50_GDSC2 <- readRDS('ic50.rds')
  library(pROC)
  ouTab <- data.frame()
  for (j in sets) {
    names(ic50_GDSC2[[j]]) <- c(j,drug_names)
    RS <- res$riskscore$`StepCox[forward] + plsRcox`[[j]]$RS
    names(RS) <- res$riskscore$`StepCox[forward] + plsRcox`[[j]]$ID
    group <- ifelse(RS>median(RS),1,0)
    ic50_GDSC2[[j]]$group <- group
    for (i in drug_names) {
      auc<-round(auc(ic50_GDSC2[[j]]$group,ic50_GDSC2[[j]][,i],direction='<'),2)
      ouTab <- rbind(
        ouTab,
        cbind(
          AUC=auc,Drug=i,Datasets=j))
    }
  }
  library(dplyr)
  library(data.table)
  tmp <- ouTab %>%
    mutate(AUC=as.numeric(AUC)) %>% 
    dplyr::filter(AUC>0.65)   #AZD1208_1449
  fwrite(tmp,file = 'ic50_selected.csv')
 
### ic50_plot ####

  library(ggpubr)
  library(stringr)
  ic50_GDSC2 <- readRDS('ic50.rds')
  drug_sel <- factor(c('AZD1208_1449'))
  my_comparisons = list( c("High", "Low"))
  for (j in sets) {
    names(ic50_GDSC2[[j]]) <- c(j,drug_names)
    RS <- res$riskscore$`StepCox[forward] + plsRcox`[[j]]$RS
    names(RS) <- res$riskscore$`StepCox[forward] + plsRcox`[[j]]$ID
    group <- factor(ifelse(RS>median(RS),'High','Low'),levels=c('Low','High'))
    pdata <- cbind(ic50_GDSC2[[j]],group)
    rocobj<-roc(pdata$group,pdata[,i],direction='<',smooth=T)
    auc<-round(auc(pdata$group,pdata[,i],direction='<',smooth=T),2)
    ggroc(rocobj,color = "red",#
            linetype = 1,
            size = 1,
            alpha = 1,
            legacy.axes = T)+#
        geom_abline(intercept = 0,#
                    slope = 1,#
                    color = "grey",
                    size = 1,
                    linetype = 1)+#
        labs(x = "False Postive Rate(1 - Specificity)",
             y = "True Positive Rate(Sensivity or Recall)")+
        annotate("text",x = 0.70,y = 0.30,#
                 label = paste("AUC =",auc),#
                 size = 8,
                 family = "serif")+#
        coord_cartesian(xlim = c(0,1),
                        ylim = c(0,1))+#
        theme_bw()+
        theme(panel.background = element_rect(fill = "transparent"),#
              axis.ticks.length = unit(0.4,"lines"),#
              axis.ticks = element_line(color = "black"),
              axis.line = element_line(size = 0.5,
                                       colour = "black"),
              axis.title = element_text(colour = "black",
                                        size = 10,
                                        face = "bold"),
              axis.text = element_text(colour = "black",
                                       size = 10,
                                       face = "bold"),
              text = element_text(size = 8,
                                  color = "black",
                                  family = "serif"))#设置文本字体
      ggsave(paste0('plot/',i,'_',j,'_AUC.pdf'),width = 5,height = 5)
    
  }
  
  