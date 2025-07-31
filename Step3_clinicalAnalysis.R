#---CGGA325----
rm(list = ls())
load("~/R/LGG/Rawdata/LGG_CGGA_325_exp_pheno.Rdata")
new_dat <- readRDS('~/R/LGG/Rawdata/CGGA_SubtypeRiskScore_Surv.rds')
pheno <- pheno[match(new_dat$sample,pheno$CGGA_ID),]
clin <- pheno %>% 
  mutate(patient=new_dat$sample,
         event=new_dat$event,
         time=new_dat$time,
         Grade=case_when(pheno$Grade=='WHO III'~1,
                          pheno$Grade=='WHO II'~0),
         Age=case_when(pheno$Age>=45~1,
                       pheno$Age<45~0),
         Gender=case_when(pheno$Gender=='Male'~1,
                         pheno$Gender=='Female'~0),
         X1p_19q=case_when(pheno$X1p19q_codeletion_status=='Non-codel'~1,
                           pheno$X1p19q_codeletion_status=='Codel'~0),
         MGMT_status=case_when(pheno$MGMTp_methylation_status=='un-methylated'~1,
                               pheno$MGMTp_methylation_status=='methylated'~0),
         IDH_status=case_when(pheno$IDH_mutation_status=='Mutant'~1,
                              pheno$IDH_mutation_status=='Wildtype'~0),
         RiskScore=new_dat$RiskScore,
         group=case_when(new_dat$RiskGroup=='High'~'high',
                         new_dat$RiskGroup=='Low'~'low'),
         Radiaotherapy=case_when(as.numeric(pheno$Radio_status..treated.1.un.treated.0.)==1~1,
                                 as.numeric(pheno$Radio_status..treated.1.un.treated.0.)==0~0)) %>% 
  select(patient,event,time,
         Gender,Age,Grade,
        X1p_19q,MGMT_status,
        IDH_status,RiskScore,
        group,Radiaotherapy)
saveRDS(clin,file = 'CGGA325_num_clin.rds')


#---CGGA693----
rm(list = ls())
load("~/R/LGG/Rawdata/CGGA_693_FPKM_LGG.Rdata")
new_dat <- readRDS('~/R/LGG/Rawdata/CGGA693_SubtypeRiskScore_Surv.rds')
pheno <- LGG_pheno[match(new_dat$sample,LGG_pheno$CGGA_ID),]
clin <- pheno %>% 
  mutate(patient=new_dat$sample,
         event=new_dat$event,
         time=new_dat$time,
         Grade=case_when(pheno$Grade=='WHO III'~1,
                         pheno$Grade=='WHO II'~0),
         Age=case_when(pheno$Age>=45~1,
                       pheno$Age<45~0),
         Gender=case_when(pheno$Gender=='Male'~1,
                          pheno$Gender=='Female'~0),
         X1p_19q=case_when(pheno$X1p19q_codeletion_status=='Non-codel'~1,
                           pheno$X1p19q_codeletion_status=='Codel'~0),
         MGMT_status=case_when(pheno$MGMTp_methylation_status=='un-methylated'~1,
                               pheno$MGMTp_methylation_status=='methylated'~0),
         IDH_status=case_when(pheno$IDH_mutation_status=='Mutant'~1,
                              pheno$IDH_mutation_status=='Wildtype'~0),
         RiskScore=new_dat$RiskScore,
         group=case_when(new_dat$RiskGroup=='High'~'high',
                         new_dat$RiskGroup=='Low'~'low'),
         Radiaotherapy=case_when(as.numeric(pheno$Radio_status..treated.1.un.treated.0.)==1~1,
                                 as.numeric(pheno$Radio_status..treated.1.un.treated.0.)==0~0)) %>% 
  select(patient,event,time,
         Gender,Age,Grade,
         X1p_19q,MGMT_status,
         IDH_status,RiskScore,
         group,Radiaotherapy)
saveRDS(clin,file = 'CGGA693_num_clin.rds')


####
library(survival)
library(data.table)
library(ggh4x)
library(dplyr)
source('clin_cox.R')
new_input_sets <- readRDS('model_input.rds')
res <- readRDS('res_ML.rds')
sets <- c('TCGA','CGGA325','CGGA301','CGGA693','GSE16011','MTAB3892') 
clin_vars <- c('Age','Gender','Grade','IDH','PQ') 

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
#vars <- c('RS','Age','Gender','Grade','IDH','PQ') 
clin <- as.data.frame(
  clin_sets[[1]]
) %>% 
  rbind(clin_sets[[2]],clin_sets[[3]],
        clin_sets[[4]],clin_sets[[5]],clin_sets[[6]])

clin$datasets <- c(rep('TCGA',498),rep('CGGA325',172),rep('CGGA301',159),
                   rep('CGGA693',420),rep('GSE16011',98),rep('MTAB3892',142))

clin <- clin %>% 
  group_by(datasets) %>% 
  mutate(group=case_when(RS>median(RS)~'High',
                   RS<= median(RS)~'Low'))
clin <- clin %>% 
  mutate(Gender=case_when(Gender=='Male'~1,
                          Gender=='Female'~0),
         Grade=case_when(Grade=='II'~1,
                         Grade=='III'~2),
         IDH=case_when(IDH=='Mut'~0,
                       IDH=='WT'~1),
         PQ=case_when(PQ=='C'~0,
                      PQ=='NC'~1),
         group=case_when(group=='High'~1,
                         group=='Low'~0))

### uni-cox #####
library(survival)
uniTab=data.frame()
for(i in colnames(clin[,4:ncol(clin)])){
  #cox <- coxph(Surv(OS.time, OS) ~ clin[,i], data = clin)
  cox <- coxph(as.formula(paste('Surv(OS.time, OS) ~ ',i )), data = clin)
  
  coxSummary = summary(cox)
  uniTab=rbind(uniTab,
               cbind(id=i,
                     HR=round(coxSummary$conf.int[,"exp(coef)"],3),
                     HR.95L=round(coxSummary$conf.int[,"lower .95"],3),
                     HR.95H=round(coxSummary$conf.int[,"upper .95"],3),
                     pvalue=round(coxSummary$coefficients[,"Pr(>|z|)"],4))
  )
}
uniTab
save(clin,file = 'LGG_clinical.Rdata')
save(uniTab,file = 'LGG_UnicoxResult.Rdata')
fwrite(uniTab,file = 'LGG_UnicoxResult.csv')
#plot
rm(list = ls())
library(forestplot)
library(dplyr)
load('LGG_UnicoxResult.Rdata')
unicox <- uniTab
unicox <- unicox[-7:-11,]
unicox <- unicox[-1,]
unicox$HR <- as.numeric(unicox$HR)
unicox$HR.95L <- as.numeric(unicox$HR.95L)
unicox$HR.95H <- as.numeric(unicox$HR.95H)
unicox$pvalue <- as.numeric(unicox$pvalue)

unicox$p.value<-ifelse(unicox$pvalue<0.001,'***',
                       ifelse(unicox$pvalue<0.01,'**',
                              ifelse(unicox$pvalue<0.05,'*','ns')))
unicox$HR2<-paste(round(unicox$HR,digits = 3)," (",
                  round(unicox$HR.95L,digits = 3),",",round(unicox$HR.95H,digits = 3),")",sep = "")
unicox<-arrange(unicox,HR)
unicox$Name <- unicox$id

pdf(file  = './Figure8A.pdf',width = 8,height = 6)
forestplot(
  unicox[,c(8,6)],
  mean = unicox[, 2],
  lower = unicox[, 3],
  upper = unicox[, 4],
  zero = 1,
  boxsize = 0.3,
  col = fpColors(box = '#CC0000', lines = 'black',
                 zero = 'grey',summary = '#99FFFF'),
  xlog=T,
  lty.ci = "dashed",
  graph.pos = 2,
  align = 'l',is.summary=F,
  hrzl_lines = list(
    "1" = gpar(lty=1),
    "7"= gpar(lty=1)),
  colgap = unit(5,'mm')
)
dev.off()
#---Figure8B-----
rm(list = ls())
library(survival)
library(survminer)
library(data.table)
load('LGG_clinical.Rdata')
clin <- clin[!duplicated(clin$ID)==T,]
clin <- clin %>% 
#  filter(!duplicated(ID)) %>% 
  tibble::column_to_rownames('ID') %>% 
  select(-datasets,-RS)

res.cox<-coxph(as.formula(paste('Surv(OS.time,OS)',sep = '~',paste(colnames(clin)[3:ncol(clin)],collapse = "+"))), 
               data = clin)
ggforest(res.cox, data = clin, 
         main = "Hazard ratio", 
         cpositions = c(0.05, 0.22, 0.4), 
         fontsize = 1, 
         refLabel = "1", noDigits = 2)
ggsave(filename = 'Figure8B.pdf',width = 8,height = 6)

#----AUC_Figure8C-E------
rm(list = ls())
library(timeROC)
load('LGG_clinical.Rdata')
load('Nomogram_model.Rdata')
clin <- clin[!duplicated(clin$ID)==T,]

clin <- clin %>% 
  #  filter(!duplicated(ID)) %>% 
  tibble::column_to_rownames('ID') %>% 
  mutate(OS.time=OS.time/365) %>% 
  select(-datasets,-RS) %>% 
  na.omit() %>% 
  as.data.frame()
clin$Nomogram <- predict(coxph(formula =  as.formula(paste("Surv(OS.time, OS) ~ ",
                                  paste(colnames(clin)[3:ncol(clin)],
                                        collapse = "+"))),
      data=clin),newdata = clin)
if(T){
  Nomogram <-with(clin, timeROC(T=OS.time,
                                delta=OS,
                                marker=Nomogram,
                                cause=1,
                                times=1,
                                iid = TRUE))
  RiskGroup <-with(clin, timeROC(T=OS.time,
                                 delta=OS,
                                 marker=group,
                                 cause=1,
                                 times=1,
                                 iid = TRUE))
  #Gender
  Gender <-with(clin, timeROC(T=OS.time,
                              delta=OS,
                              marker=Gender,
                              cause=1,
                              times=1,
                              iid = TRUE))
  #Age
  Age <-with(clin, timeROC(T=OS.time,
                           delta=OS,
                           marker=Age,
                           cause=1,
                           times=1,
                           iid = TRUE))
  #Grade
  Grade <-with(clin, timeROC(T=OS.time,
                             delta=OS,
                             marker=Grade,
                             cause=1,
                             times=1,
                             iid = TRUE))
  #X1p_19q
  X1p_19q <-with(clin, timeROC(T=OS.time,
                               delta=OS,
                               marker=PQ,
                               cause=1,
                               times=1,
                               iid = TRUE))
  
  #IDH_status
  IDH_status <-with(clin, timeROC(T=OS.time,
                                  delta=OS,
                                  marker=IDH,
                                  cause=1,
                                  times=1,
                                  iid = TRUE))
  
  pdf(file = 'Figure8C.pdf',width = 6,height = 6)
  plot(Nomogram, time = 1, col="#92C5DE", lwd=2, title = "")
  plot(RiskGroup, time = 1, col="#E41A1C", lwd=2, add = T)
  plot(Gender, time = 1, col="grey", lwd=2, add = T)
  plot(Age, time = 1, col="grey", lwd=2, add = T)
  plot(Grade, time = 1, col="grey", lwd=2, add = T)
  plot(X1p_19q, time = 1, col="grey", lwd=2, add = T)
  plot(IDH_status, time = 1, col="grey", lwd=2, add = T)
  legend("bottomright",
         c(paste0("Nomogram: ",round(Nomogram[["AUC"]][2],2)), 
           paste0("RiskGroup: ",round(RiskGroup[["AUC"]][2],2)), 
           paste0("Gender: ",round(Gender[["AUC"]][2],2)), 
           paste0("Age: ",round(Age[["AUC"]][2],2)),
           paste0("Grade: ",round(Grade[["AUC"]][2],2)),
           paste0("X1p_19q: ",round(X1p_19q[["AUC"]][2],2)),
           paste0("IDH_status: ",round(IDH_status[["AUC"]][2],2))
           
         ),
         col=c("#92C5DE","#E41A1C", "grey", "grey", "grey", 
               "grey", "grey"),
         lty=1, lwd=2,bty = "n",title.font =30
  )
  title('1-years of AUC')
  dev.off()
}
#3year'
if(T){
  Nomogram <-with(clin, timeROC(T=OS.time,
                                delta=OS,
                                marker=Nomogram,
                                cause=1,
                                times=3,
                                iid = TRUE))
  RiskGroup <-with(clin, timeROC(T=OS.time,
                                 delta=OS,
                                 marker=group,
                                 cause=1,
                                 times=3,
                                 iid = TRUE))
  #Gender
  Gender <-with(clin, timeROC(T=OS.time,
                              delta=OS,
                              marker=Gender,
                              cause=1,
                              times=3,
                              iid = TRUE))
  #Age
  Age <-with(clin, timeROC(T=OS.time,
                           delta=OS,
                           marker=Age,
                           cause=1,
                           times=3,
                           iid = TRUE))
  #Grade
  Grade <-with(clin, timeROC(T=OS.time,
                             delta=OS,
                             marker=Grade,
                             cause=1,
                             times=3,
                             iid = TRUE))
  #X1p_19q
  X1p_19q <-with(clin, timeROC(T=OS.time,
                               delta=OS,
                               marker=PQ,
                               cause=1,
                               times=3,
                               iid = TRUE))
  #IDH_status
  IDH_status <-with(clin, timeROC(T=OS.time,
                                  delta=OS,
                                  marker=IDH,
                                  cause=1,
                                  times=3,
                                  iid = TRUE))
  
  pdf(file = 'Figure8D.pdf',width = 6,height = 6)
  plot(Nomogram, time = 3, col="#92C5DE", lwd=2, title = "")
  plot(RiskGroup, time = 3, col="#E41A1C", lwd=2, add = T)
  plot(Gender, time = 3, col="grey", lwd=2, add = T)
  plot(Age, time = 3, col="grey", lwd=2, add = T)
  plot(Grade, time = 3, col="grey", lwd=2, add = T)
  plot(X1p_19q, time = 3, col="grey", lwd=2, add = T)
  plot(IDH_status, time = 3, col="grey", lwd=2, add = T)
  legend("bottomright",
         c(paste0("Nomogram: ",round(Nomogram[["AUC"]][2],2)), 
           paste0("RiskGroup: ",round(RiskGroup[["AUC"]][2],2)), 
           paste0("Gender: ",round(Gender[["AUC"]][2],2)), 
           paste0("Age: ",round(Age[["AUC"]][2],2)),
           paste0("Grade: ",round(Grade[["AUC"]][2],2)),
           paste0("X1p_19q: ",round(X1p_19q[["AUC"]][2],2)),
           paste0("IDH_status: ",round(IDH_status[["AUC"]][2],2))
         ),
         col=c("#92C5DE","#E41A1C", "grey", "grey", "grey", 
               "grey", "grey"),
         lty=1, lwd=2,bty = "n"
  )
  title('3-years of AUC')
  dev.off()
}
#5year
if(T){
  Nomogram <-with(clin, timeROC(T=OS.time,
                                delta=OS,
                                marker=Nomogram,
                                cause=1,
                                times=5,
                                iid = TRUE))
  RiskGroup <-with(clin, timeROC(T=OS.time,
                                 delta=OS,
                                 marker=group,
                                 cause=1,
                                 times=5,
                                 iid = TRUE))
  #Gender
  Gender <-with(clin, timeROC(T=OS.time,
                              delta=OS,
                              marker=Gender,
                              cause=1,
                              times=5,
                              iid = TRUE))
  #Age
  Age <-with(clin, timeROC(T=OS.time,
                           delta=OS,
                           marker=Age,
                           cause=1,
                           times=5,
                           iid = TRUE))
  #Grade
  Grade <-with(clin, timeROC(T=OS.time,
                             delta=OS,
                             marker=Grade,
                             cause=1,
                             times=5,
                             iid = TRUE))
  #X1p_19q
  X1p_19q <-with(clin, timeROC(T=OS.time,
                               delta=OS,
                               marker=PQ,
                               cause=1,
                               times=5,
                               iid = TRUE))
  
  #IDH_status
  IDH_status <-with(clin, timeROC(T=OS.time,
                                  delta=OS,
                                  marker=IDH,
                                  cause=1,
                                  times=5,
                                  iid = TRUE))
  
  
  pdf(file = 'Figure8E.pdf',width = 6,height = 6)
  plot(Nomogram, time = 5, col="#92C5DE", lwd=2, title = "")
  plot(RiskGroup, time = 5, col="#E41A1C", lwd=2, add = T)
  plot(Gender, time = 5, col="grey", lwd=2, add = T)
  plot(Age, time = 5, col="grey", lwd=2, add = T)
  plot(Grade, time = 5, col="grey", lwd=2, add = T)
  plot(X1p_19q, time = 5, col="grey", lwd=2, add = T)
  plot(IDH_status, time = 5, col="grey", lwd=2, add = T)
  legend("bottomright",
         c(paste0("Nomogram: ",round(Nomogram[["AUC"]][2],2)), 
           paste0("RiskGroup: ",round(RiskGroup[["AUC"]][2],2)), 
           paste0("Gender: ",round(Gender[["AUC"]][2],2)), 
           paste0("Age: ",round(Age[["AUC"]][2],2)),
           paste0("Grade: ",round(Grade[["AUC"]][2],2)),
           paste0("X1p_19q: ",round(X1p_19q[["AUC"]][2],2)),
           paste0("IDH_status: ",round(IDH_status[["AUC"]][2],2))
         ),
         col=c("#92C5DE","#E41A1C", "grey", "grey", "grey", 
               "grey", "grey"),
         lty=1, lwd=2,bty = "n",title.font =30
  )
  title('5-years of AUC')
  dev.off()
}

#---Nomogram_Figure8F-----
setwd('~/R/LGG/')
rm(list = ls())
library(glmnet)
library(rms)
library(dplyr)
load('LGG_clinical.Rdata')
clin <- clin[!duplicated(clin$ID)==T,]
clin <- clin %>% 
  tibble::column_to_rownames('ID') %>% 
  select(-datasets,-RS)

clin <- clin %>% 
  mutate(Grade=factor(case_when(Grade==2~'WHO III',
                         Grade==1~'WHO II')),
         Age=factor(case_when(Age==1~'>40',
                       Age==0~'<40')),
         Gender=factor(case_when(Gender==1~'Male',
                          Gender==0~'Female')),
         X1p_19q=factor(case_when(PQ==1~'Non-codel',
                           PQ==0~'Codel')),
         IDH=factor(case_when(IDH==0~'Mutant',
                              IDH==1~'Wildtype')),
         group=factor(case_when(group==1~'High',
                         group==0~'Low'))) %>% 
  select(OS.time,OS,
         Age,Grade,
         X1p_19q,IDH,group)
saveRDS(clin,file = 'LGG_clinical_summary.rds')

#Nomogram_Plot----
rm(list = ls())
library(rms)
load('LGG_clinical.Rdata')
clin <- clin[!duplicated(clin$ID)==T,]
clin <- clin %>% 
  #  filter(!duplicated(ID)) %>% 
  tibble::column_to_rownames('ID') %>% 
  mutate(OS.time=OS.time/365) %>% 
  select(-datasets,-RS,-Gender)

clin <- readRDS('LGG_clinical_summary.rds') %>% 
  mutate(OS.time=OS.time/365) 
ph <- na.omit(clin) %>% 
  as.data.frame()
dd<-datadist(ph)
options(datadist="dd")
mod <- cph(formula =  as.formula(paste("Surv(OS.time, OS) ~ ",
                                       paste(colnames(ph)[3:ncol(ph)],
                                             collapse = "+"))),
           data=ph,x=T,y=T,surv = T)

surv<-Survival(mod) 
surv1<-function(x) surv(1,x)
surv3<-function(x) surv(3,x)
surv5<-function(x) surv(5,x)

pdf(file = './Figure8F.pdf',width = 10,height = 8)
x<-nomogram(mod,
            fun = list(surv1,surv3,surv5),
            funlabel = c('1-year survival Probability',
                         '3-year survival Probability',
                         '5-year survival Probability'))

plot(x)
dev.off()
save(mod,file = 'Nomogram_model.Rdata')
#----Figure8G----
rm(list=ls())
library(dplyr)
library(survival)
library(rms) #2023 v6.7-1

library(pec)
clin <- readRDS('LGG_clinical_summary.rds') %>% 
  mutate(OS.time=OS.time/365)

rt <- clin 
rt <- na.omit(rt)
bioCol=rainbow(ncol(rt)-2, s=0.9, v=0.9)

dd<-datadist(rt)
options(datadist="dd")
#plot
Nomogram <- cph(formula =  as.formula(paste("Surv(OS.time, OS) ~ ",paste(colnames(rt)[3:ncol(rt)],collapse = "+"))),
                data=rt,x=T,y=T,surv = T)
#Gender <- cph(Surv(time, event)~Gender, data=rt, surv=TRUE)
Age <- cph(Surv(OS.time, OS)~Age, data=rt, surv=TRUE)
Grade <- cph(Surv(OS.time, OS)~Grade, data=rt, surv=TRUE)
X1p_19q <- cph(Surv(OS.time, OS)~X1p_19q, data=rt, surv=TRUE)
IDH_status <- cph(Surv(OS.time, OS)~IDH, data=rt, surv=TRUE)
RiskGroup=cph(Surv(OS.time, OS)~group, data=rt, surv=TRUE)
c_index  <- cindex(list('Nomogram'=Nomogram,
                        "Group"=RiskGroup,
                        'Age'=Age,
                        '1p/19q codeleted'=X1p_19q,
                        'IDH status'=IDH_status,
                        'Grade'=Grade),
                   formula=Surv(OS.time, OS)~ .,data=rt,eval.times=seq(0,50,1),
                   splitMethod="bootcv",B=1000)
pdf(file="Figure8G.pdf", width=5.5, height=5)
plot(c_index, 
     xlim=c(1,10), ylim=c(0.4,0.8), 
     col=bioCol, xlab="Time (years)",
     legend.x=8, legend.y=0.93, legend.cex=0.8)
dev.off()
#--Figure8H-----
rm(list=ls())
library(dplyr)
library(glmnet)
library(rms)
load('LGG_clinical.Rdata')
clin <- clin[!duplicated(clin$ID)==T,]
clin <- clin %>% 
 
  tibble::column_to_rownames('ID') %>% 
  mutate(OS.time=OS.time/365) %>% 
  select(-datasets,-RS,-Gender)
ph <- na.omit(clin)
dd<-datadist(ph) 
options(datadist="dd") #
if(T){
  f1 <- cph(formula =  as.formula(paste("Surv(OS.time, OS) ~ ",paste(colnames(ph)[3:ncol(ph)],collapse = "+"))),
            data=ph,x=T,y=T,surv = T, time.inc=1)
  cal1 <- calibrate(f1, cmethod="KM", method="boot", u=1, m=260, B=950)
  f3 <- cph(formula =  as.formula(paste("Surv(OS.time, OS) ~ ",paste(colnames(ph)[3:ncol(ph)],collapse = "+"))),
            data=ph,x=T,y=T,surv = T, time.inc=3)
  cal3 <- calibrate(f3, cmethod="KM", method="boot", u=3, m=260, B=950)
  
  f5 <- cph(formula =  as.formula(paste("Surv(OS.time, OS) ~ ",paste(colnames(ph)[3:ncol(ph)],collapse = "+"))),
            data=ph,x=T,y=T,surv = T,  time.inc=5)
  cal5 <- calibrate(f5, cmethod="KM", method="boot", u=5, m=260, B=950)
  pdf(file = 'Figure8H.pdf',width = 6,height = 6)
  
  plot(cal1,lwd = 2,lty = 0,errbar.col = c("#7FFFD4"),
       bty = "l", 
       xlim = c(0,1),ylim= c(0,1),
       xlab = "Nomogram-predicted OS (%)",ylab = "Observed OS (%)",
       col = c("#7FFFD4"),
       cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6)
  lines(cal1[,c('mean.predicted',"KM")],
        type = 'b', lwd = 1, col = c("#7FFFD4"), pch = 16)
  mtext("")
  plot(cal3,lwd = 2,lty = 0,errbar.col = c("#2166AC"),
       xlim = c(0,1),ylim= c(0,1),
       col = c("#2166AC"),
       cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6,add = T)
  lines(cal3[,c('mean.predicted',"KM")],
        type = 'b', lwd = 1, col = c("#2166AC"), pch = 16)
  mtext("")
  plot(cal5,lwd = 2,lty = 0,errbar.col = c("#B2182B"),
       xlim = c(0,1),ylim= c(0,1),col = c("#B2182B"),add = T)
  lines(cal5[,c('mean.predicted',"KM")],
        type = 'b', lwd = 1, col = c("#B2182B"), pch = 16)
  abline(0,1, lwd = 2, lty = 3, col = c("#224444"))
  
  legend("topleft", 
         legend = c("1-year","3-year","5-year"),
         col =c("#7FFFD4","#2166AC","#B2182B"), 
         lwd = 2,
         cex = 1.2,
         bty = "n")
  dev.off()
}
