rm(list=ls())
library(tinyarray)
library(readxl)
library(clusterProfiler)
library(data.table)
library(stringr)
library(dplyr)

#### TCGA-LGG #####
rm(list=ls())
library(IOBR)
library(data.table)

load('D:\\SynologyDrive\\R\\data\\TCGA-LGG_conut_expr_clin.Rdata')
exp <- exp[,substr(colnames(exp),15,16)=='1A']
colnames(exp) <- substr(colnames(exp),1,12)
surv <- surv[match(colnames(exp),surv$`_PATIENT`,nomatch = F),]
surv <- surv %>% 
  filter(is.na(OS.time)!=T,is.na(OS)!=T)
exp <- exp[,match(surv$`_PATIENT`,colnames(exp),nomatch = F)]
clinical <- clinical[match(colnames(exp),clinical$bcr_patient_barcode,nomatch = F),]
rownames(clinical) <- NULL
LGG_clinical_summary <- LGG_clinical_summary[match(colnames(exp),rownames(LGG_clinical_summary),nomatch = T),]

clin <- clinical %>% 
  column_to_rownames(var = 'bcr_patient_barcode') %>% 
  mutate(OS.time=as.double(surv$OS.time),
         OS=as.double(surv$OS),
         Gender=case_when(gender=='FEMALE'~'Female',
                          gender=='MALE'~'Male'),
         PQ=case_when(LGG_clinical_summary$X1p_19q=='Non-codel'~'NC',
                      LGG_clinical_summary$X1p_19q=='Codel'~'C'),
         MGMT=case_when(LGG_clinical_summary$MGMT_status=='un-methylated'~'NM',
                        LGG_clinical_summary$MGMT_status=='methylated'~'M'),
         Age=as.numeric(age_at_initial_pathologic_diagnosis),
         IDH=case_when(LGG_clinical_summary$IDH_status=='Wildtype'~'WT',
                       LGG_clinical_summary$IDH_status=='Mutant'~'Mut'),
         Grade=case_when(LGG_clinical_summary$Grade=='WHO III'~'III',
                         LGG_clinical_summary$Grade=='WHO II'~'II'),
         Radio=LGG_clinical_summary$Radiotherapy,
         Histology=case_when(histological_type=='Astrocytoma'~'A',
                             histological_type=='Oligoastrocytoma'~'OA',
                             histological_type=='Oligodendroglioma'~'OD')) %>% 
  select(OS.time,OS,Age,Gender,Grade,Histology,IDH,PQ,MGMT,Radio)
saveRDS(exp,file = 'cleandata/TCGA_exp.rds')
saveRDS(clin,file = 'cleandata/TCGA_clin.rds')
save(exp,clin,file = 'TCGA_exp_clin.Rdata')
exp <- readRDS('cleandata/TCGA_exp.rds')
clin <- readRDS('cleandata/TCGA_clin.rds')

exp <- count2tpm(exp,idType='Symbol',gene_symbol = "symbol")
#### CGGA_325 #####
rm(list = ls())
exp <- fread('data/CGGA.mRNAseq_325.Read_Counts-genes.20220620.txt',
             stringsAsFactors = F,data.table = F) %>% 
  tibble::column_to_rownames('gene_name')
clinical <- fread('data/CGGA.mRNAseq_325_clinical.20200506.txt')
str(clinical)
clin <- clinical %>% 
  tibble::column_to_rownames('CGGA_ID') %>% 
  mutate(OS.time=as.double(OS)) %>% 
  mutate(OS=as.double(`Censor (alive=0; dead=1)`),
         IDH=case_when(IDH_mutation_status=='Wildtype'~'WT',
                       IDH_mutation_status=='Mutant'~'Mut'),
         Grade=case_when(Grade=='WHO IV'~'IV',
                         Grade=='WHO III'~'III',
                         Grade=='WHO II'~'II'),
         PQ=case_when(`1p19q_codeletion_status`=='Non-codel'~'NC',
                      `1p19q_codeletion_status`=='Codel'~'C'),
         MGMT=case_when(`MGMTp_methylation_status`=='un-methylated'~'NM',
                      `MGMTp_methylation_status`=='methylated'~'M')) %>% 
  dplyr::rename(Radio=`Radio_status (treated=1;un-treated=0)`,
                Chemo=`Chemo_status (TMZ treated=1;un-treated=0)`) %>% 
  filter(Grade=='II'|Grade=='III') %>% 
  filter(is.na(OS.time)!=T) %>% 
  dplyr::select(OS.time,OS,Age,Gender,Grade,Histology,IDH,PQ,MGMT,Radio,Chemo,PRS_type)
exp <- exp[,match(rownames(clin),colnames(exp))]
exp <- count2tpm(exp,idType='Symbol')
saveRDS(exp,file = 'cleandata/CGGA325_exp.rds')
saveRDS(clin,file = 'cleandata/CGGA325_clin.rds')
save(exp,clin,file = 'CGGA325_exp_clin.Rdata')
#### CGGA_693 #####
rm(list = ls())
exp <- fread('data/CGGA.mRNAseq_693.Read_Counts-genes.20220620.txt',
             stringsAsFactors = F,data.table = F) %>% 
  tibble::column_to_rownames('gene_name')
clinical <- fread('data/CGGA.mRNAseq_693_clinical.20200506.txt')
str(clinical)
clin <- clinical %>% 
  tibble::column_to_rownames('CGGA_ID') %>% 
  mutate(OS.time=as.double(OS)) %>% 
  mutate(OS=as.double(`Censor (alive=0; dead=1)`),
         IDH=case_when(IDH_mutation_status=='Wildtype'~'WT',
                       IDH_mutation_status=='Mutant'~'Mut'),
         Grade=case_when(Grade=='WHO IV'~'IV',
                         Grade=='WHO III'~'III',
                         Grade=='WHO II'~'II'),
         Histology=case_when(Histology=='A'~'A',
                             Histology=='rA'~'A',
                             Histology=='AA'~'AA',
                             Histology=='rAA'~'AA',
                             Histology=='AO'~'AO',
                             Histology=='rAO'~'AO',
                             Histology=='AOA'~'rAOA',
                             Histology=='AOA'~'rAOA',
                             Histology=='O'~'rO',
                             Histology=='O'~'rO',
                             Histology=='OA'~'rOA',
                             Histology=='OA'~'rOA'),
         PQ=case_when(`1p19q_codeletion_status`=='Non-codel'~'NC',
                      `1p19q_codeletion_status`=='Codel'~'C'),
         MGMT=case_when(`MGMTp_methylation_status`=='un-methylated'~'NM',
                        `MGMTp_methylation_status`=='methylated'~'M')) %>% 
  dplyr::rename(Radio=`Radio_status (treated=1;un-treated=0)`,
                Chemo=`Chemo_status (TMZ treated=1;un-treated=0)`) %>% 
  filter(Grade=='II'|Grade=='III') %>% 
  filter(is.na(OS.time)!=T) %>% 
  dplyr::select(OS.time,OS,Age,Gender,Grade,Histology,IDH,PQ,MGMT,Radio,Chemo,PRS_type)
exp <- exp[,match(rownames(clin),colnames(exp))]
exp <- count2tpm(exp,idType='Symbol')
saveRDS(exp,file = 'cleandata/CGGA693_exp.rds')
saveRDS(clin,file = 'cleandata/CGGA693_clin.rds')
save(exp,clin,file = 'CGGA693_exp_clin.Rdata')

#### CGGA_301 #####
rm(list = ls())
exp <- fread('data/CGGA.mRNA_array_301_gene_level.20200506.txt',
             stringsAsFactors = F,data.table = F) %>% 
  tibble::column_to_rownames('Gene_Name')
clinical <- fread('data/CGGA.mRNA_array_301_clinical.20200506.txt')
str(clinical)
clin <- clinical %>% 
  tibble::column_to_rownames('CGGA_ID') %>% 
  mutate(OS.time=as.double(OS)) %>% 
  mutate(OS=as.double(`Censor (alive=0; dead=1)`),
         IDH=case_when(IDH_mutation_status=='Wildtype'~'WT',
                       IDH_mutation_status=='Mutant'~'Mut'),
         Grade=case_when(Grade=='WHO IV'~'IV',
                         Grade=='WHO III'~'III',
                         Grade=='WHO II'~'II'),
         PQ=case_when(`1p19q_Codeletion_status`=='Non-codel'~'NC',
                      `1p19q_Codeletion_status`=='Codel'~'C'),
         MGMT=case_when(`MGMTp_methylation_status`=='un-methylated'~'NM',
                        `MGMTp_methylation_status`=='methylated'~'M')) %>% 
  dplyr::rename(Radio=`Radio_status (treated=1;un-treated=0)`,
                Chemo=`Chemo_status (TMZ treated=1;un-treated=0)`) %>% 
  mutate(Type=case_when(Grade=='IV'~'GBM',
                        Grade=='III'~'LGG',
                        Grade=='II'~'LGG')) %>% 
  filter(Type=='LGG') %>% 
  filter(is.na(OS.time)!=T) %>% 
  dplyr::select(OS.time,OS,Age,Gender,Grade,Histology,IDH,PQ,MGMT,Radio,Chemo,PRS_type)
exp <- exp[,match(rownames(clin),colnames(exp))]

saveRDS(exp,file = 'cleandata/CGGA301_exp.rds')
saveRDS(clin,file = 'cleandata/CGGA301_clin.rds')
save(exp,clin,file = 'CGGA301_exp_clin.Rdata')

#### GSE108474 #####
rm(list = ls())
library(tinyarray)
library(GEOquery)
gset <- getGEO(GEO = 'GSE108474',AnnotGPL = F,getGPL = F)
ids <- AnnoProbe::idmap('GPL570')
dat <- gset$GSE108474_series_matrix.txt.gz@assayData$exprs
exp <- trans_array(dat,ids) %>% 
  as.data.frame()
phe <- gset$GSE108474_series_matrix.txt.gz@phenoData@data %>% 
  as.data.frame() %>% 
  dplyr::select(title,`extract name:ch1`)
clinical <- fread('data/GSE108474_REMBRANDT_clinical.data.txt') %>% 
  dplyr::filter(WHO_GRADE=='II'|WHO_GRADE=='III') %>% 
  dplyr::filter(is.na(EVENT_OS)==F) %>% 
  dplyr::filter(GENDER=='MALE'|GENDER=='FEMALE')
clinical <- clinical %>% 
  mutate(OS.time=as.double(OVERALL_SURVIVAL_MONTHS)) %>% 
  mutate(OS=as.double(EVENT_OS),
         Age=as.numeric(stringr::str_split_fixed(AGE_RANGE,pattern = '[-]',n=2)[,1]),
         Grade=WHO_GRADE,
         Histology=case_when(DISEASE_TYPE=='ASTROCYTOMA'~'A',
                             DISEASE_TYPE=='OLIGODENDROGLIOMA'~'OD'),
         Gender=case_when(GENDER=='MALE'~'Male',
                          GENDER=='FEMALE'~'Female'),
         title=SUBJECT_ID) %>% 
  dplyr::select(title,OS.time,OS,Age,Gender,Grade,Histology)

phe <- phe[match(clinical$title,phe$title,nomatch = F),]
exp <- exp[,match(rownames(phe),colnames(exp))]
clinical <- clinical[match(phe$title,clinical$title,nomatch = F),]
rownames(clinical) <- NULL
clin <- as.data.frame(clinical)
rownames(clin) <- colnames(exp)
clin <- clin[,-1]
clin$OS.time <- clin$OS.time*30
clin <- clin[is.na(clin$OS.time)!=T,]
saveRDS(exp,file = 'cleandata/GSE108474_exp.rds')
saveRDS(clin,file = 'cleandata/GSE108474_clin.rds')
save(exp,clin,file = 'GSE108474_exp_clin.Rdata')
#### GSE10611 ####
rm(list = ls())
library(tinyarray)
library(GEOquery)
library(dplyr)
library(GroundWork)
# devtools::install_github("grswsci/GroundWork")
exp <- read.table('GSE16011_series_matrix.txt',
                  header = T,comment.char = '!') %>% 
  tibble::column_to_rownames('ID_REF')
ids <- getGPL(GPL_ID = 'GPL8542')
colnames(ids) <- c('probe_id','symbol')
exp <- trans_array(exp,ids) %>% 
  as.data.frame()
Clinical <- getClinical(GEO_ID = 'GSE16011') 
colnames(Clinical) <- Clinical[1,]
Clinical <- Clinical[-1,]
Clinical <- as.data.frame(Clinical)
Clinical$Number <- as.numeric(stringr::str_split_fixed(rownames(Clinical),pattern = '[.]',2)[,2])
phe <- readxl::read_xlsx('data/GSE16011_clin.xlsx')
phe <- phe %>% 
  filter(`Reviewed histological diagnosis`!='control',is.na(event)!=T)
Clinical <- Clinical[match(phe$Number,Clinical$Number,nomatch = F),]
phe$GEO_ID <- Clinical$`!Sample_geo_accession`

clinical <- phe %>% 
  tibble::column_to_rownames('GEO_ID') %>% 
  mutate(OS=event,
         OS.time=time,
         Age=as.numeric(stringr::str_split_fixed(`Age at diagnosis`,pattern = '_',n=2)[,1]),
         Type=case_when(Grade=='IV'~'GBM',
                        Grade=='III'~'LGG',
                        Grade=='II'~'LGG')) %>% 
  filter(Type=='LGG') %>% 
  dplyr::select(OS,OS.time,Age,Gender,Grade,Histology,IDH,PQ,Radio,Chemo)
exp <- exp[,match(rownames(clinical),colnames(exp))]
clinical$OS.time <- clinical$OS.time*365
clinical <- clinical[is.na(clinical$OS.time)!=T,]

saveRDS(exp,file = 'cleandata/GSE16011_exp.rds')
saveRDS(clinical,file = 'cleandata/GSE16011_clin.rds')
save(exp,clinical,file = 'GSE16011_exp_clin.Rdata')



#### E-MTAB-3892 #####
setwd('D:\\3892\\')
library(ArrayExpress)
library(oligo)
library(dplyr)
library(hgu133plus2.db)
library(tinyarray)
#获取expr、pheno和annotation
raw_data <- read.celfiles(filenames = file.path('raw',SDRF$Array.Data.File),
                          verbose = T)
rownames(SDRF) <- sampleNames(raw_data)

SDRF <- AnnotatedDataFrame(SDRF)
phenoData(raw_data)<- SDRF
pd <- pData(raw_data)

raw_data %>%  rma(.,normalize =TRUE, subset = NULL) %>% 
  exprs(.) -> expr_raw
expr_raw <- as.data.frame(expr_raw)
SYMBOL <- annotate::lookUp(rownames(expr_raw),"hgu133plus2.db", "SYMBOL")
length(SYMBOL)
ids = data.frame(probe_id=rownames(expr_raw),symbol=as.vector(unlist(SYMBOL)),stringsAsFactors = F)
ids <- ids[is.na(ids$symbol)!=T,]
exp <- trans_array(expr_raw,ids,from = 'probe_id',to = 'symbol')
#save(exp,pd,file = 'raw_E-MTAB-3892.Rdata')
table(pd$Characteristics.histology.)
#数据筛选
clin <- pd %>% 
  filter(Characteristics.disease.staging.=='tumor',
         Factor.Value.histology.!='normal',
         Factor.Value.histology.!='NOS') %>% 
  mutate(OS.time=as.numeric(Characteristics.os.delay.)*30,
         OS=as.numeric(Factor.Value.os.event.),
         Age=as.numeric(Characteristics.age.),
         Gender=case_when(Characteristics.sex.=='male'~'Male',
                          Characteristics.sex.=='female'~'Female'),
         ID=Source.Name,
         ArrayID=Array.Data.File,
         Grade=case_when(Factor.Value.histology_grade.=='2'~'II',
                         Factor.Value.histology_grade.=='3'~'III',
                         Factor.Value.histology_grade.=='4'~'IV'),
         IDH=case_when(Factor.Value.idh1.gene.mutation.=='M'~'Mut',
                       Factor.Value.idh1.gene.mutation.=='WT'~'WT'),
         PQ=case_when(Factor.Value.1p.19q.co.deletion.=='1'~'C',
                      Factor.Value.1p.19q.co.deletion.=='0'~'NC'),
         Histology=case_when(Characteristics.histology.=='Low-grade Oligodendroglioma'~'OD',
                             Characteristics.histology.=='Low-grade Oligoastrocytoma'~'OA',
                             Characteristics.histology.=='High-grade Oligodendroglioma'~'OD',
                             Characteristics.histology.=='High-grade Oligoastrocytoma'~'OA',
                             Characteristics.histology.=='Diffuse astrocytoma'~'A')) %>% 
  dplyr::select(ID,OS.time,OS,Age,Gender,Grade,Histology,IDH,PQ,ArrayID) %>% 
  filter(is.na(OS.time)!=T,is.na(OS)!=T) %>% 
  filter(Grade=='II'|Grade=='III')
exp <- exp[,match(clin$ArrayID,colnames(exp),nomatch = F)]
colnames(exp) <- clin$ID
clin <- clin[,-10]
saveRDS(exp,file = 'D:/SynologyDrive/R/LY/LGG/cleandata/MTAB3892_exp.rds')
saveRDS(clin,file = 'D:/SynologyDrive/R/LY/LGG/cleandata/MTAB3892_clin.rds')
save(exp,clin,file = 'D:/SynologyDrive/R/LY/LGG/MTAB3892_exp_clin.Rdata')


