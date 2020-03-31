library(reshape2)
library(dplyr)
library(purrr)
library(broom)
library(survMisc)
library(survminer)
library(dplyr)
library(ggrepel)
library('XML')
library('methods')
library(ggplot2)
library(gridExtra)
library(colorRamps)
library(stringr)
library(cowplot)


#### Format clinical data ####
path = '/Volumes/BioDiscovery/Kristin/SurVarCN_previousVersions/SurVarCN_updatedTCGA/cBioPortal'
temp = list.files(path = path,pattern = '*clinical_data_patient.txt',full.names = TRUE)
clin = lapply(temp, read.delim, stringsAsFactors = F)
names(clin) = gsub('_.*','',basename(temp))

for (i in 1:length(clin)) {
  #restructure df
  colnames(clin[[i]]) = clin[[i]][4,]
  clin[[i]] = clin[[i]][-c(1:4),]
  # add missing columns and populate with NA
  req_cols = c('PATIENT_ID','HISTOLOGICAL_SUBTYPE','AGE','SEX','RACE','AJCC_PATHOLOGIC_TUMOR_STAGE','HISTOLOGICAL_DIAGNOSIS','GRADE','TREATMENT_OUTCOME_FIRST_COURSE','RADIATION_TREATMENT_ADJUVANT','OS_STATUS','OS_MONTHS','DFS_STATUS','DFS_MONTHS')
  if (!all(req_cols%in%colnames(clin[[i]]))) {
    miss_cols = req_cols[!req_cols %in% colnames(clin[[i]])]
    temp_df = data.frame(matrix(ncol = length(miss_cols),nrow = nrow(clin[[i]])))
    colnames(temp_df) = miss_cols
    clin[[i]] = cbind(clin[[i]],temp_df)
  }
  clin[[i]] = subset(clin[[i]],select=req_cols)
  colnames(clin[[i]])[colnames(clin[[i]])=='HISTOLOGICAL_SUBTYPE'] = 'SUBTYPE'
  clin[[i]][clin[[i]]=='[Not Available]'] = NA
  clin[[i]][clin[[i]]=='[Not Applicable]'] = NA
  # convert survival time from months to days
  colnames(clin[[i]])[colnames(clin[[i]])=='OS_MONTHS'] = 'OS_DAYS'
  colnames(clin[[i]])[colnames(clin[[i]])=='DFS_MONTHS'] = 'DFS_DAYS'
  clin[[i]]$OS_DAYS = round(as.numeric(clin[[i]]$OS_DAYS) * 30.42)
  clin[[i]]$DFS_DAYS = round(as.numeric(clin[[i]]$DFS_DAYS) * 30.42)
  # create binary 1/0 classifications for survival
  clin[[i]]$OS_STATUS = ifelse(clin[[i]]$OS_STATUS=='LIVING',0,1)
  clin[[i]]$DFS_STATUS = ifelse(clin[[i]]$DFS_STATUS=='DiseaseFree',0,1)
  # edit age to match TRGAted files, age range instead of exact number
  clin[[i]]$AGE = as.numeric(clin[[i]]$AGE)
  clin[[i]]$AGE2[clin[[i]]$AGE<=20]='0-20'
  clin[[i]]$AGE2[clin[[i]]$AGE>20 & clin[[i]]$AGE<=30]='21-30'
  clin[[i]]$AGE2[clin[[i]]$AGE>30 & clin[[i]]$AGE<=40]='31-40'
  clin[[i]]$AGE2[clin[[i]]$AGE>40 & clin[[i]]$AGE<=50]='41-50'
  clin[[i]]$AGE2[clin[[i]]$AGE>50 & clin[[i]]$AGE<=60]='51-60'
  clin[[i]]$AGE2[clin[[i]]$AGE>60 & clin[[i]]$AGE<=70]='61-70'
  clin[[i]]$AGE2[clin[[i]]$AGE>70 & clin[[i]]$AGE<=80]='71-80'
  clin[[i]]$AGE2[clin[[i]]$AGE>80 & clin[[i]]$AGE<=90]='81-90'
  clin[[i]]$AGE2[clin[[i]]$AGE>90 & clin[[i]]$AGE<=100]='91-100'
  clin[[i]]$AGE2[clin[[i]]$AGE>100]='>100'
  clin[[i]]$AGE = clin[[i]]$AGE2
  clin[[i]] = subset(clin[[i]],select=-AGE2)
}


#### merge subtype info into clinical data ####

# read in subtype csvs
path = '/Volumes/BioDiscovery/Kristin/SurVarCN_previousVersions/SurVarCN_updatedTCGA/SUBTYPE_downloaded'
temp = list.files(path = path,pattern = '*.csv',full.names = TRUE)
SUBs = lapply(temp, read.delim, stringsAsFactors = F, sep=',')
names(SUBs) = gsub('_subtype.*csv','',basename(temp))

# initialize new list
clinSUB = list()


### GBM ####
dat = clin[['gbm']]
# subtypes
sub = SUBs[['gbm_lgg']]
colnames(sub) = sub[1,]
sub = sub[-1,]
sub = sub[sub$Study=='Glioblastoma multiforme',]
sub = subset(sub,select=c('Case','Original Subtype'))
# merge
dat = merge(dat,sub,by.x='PATIENT_ID',by.y='Case',all.x=TRUE)
dat$SUBTYPE = dat$`Original Subtype`
dat = dat[,-ncol(dat)]
clinSUB[['gbm']] = dat


#### BRCA ####
dat = clin[['brca']]
# subtypes
sub = SUBs[['brca_cesc_ov_ucec_ucs']]
colnames(sub) = sub[1,]
sub = sub[-1,]
sub = subset(sub,select=c(Sample.ID,BRCA_Subtype_PAM50))
# merge
dat = merge(dat,sub,by.x='PATIENT_ID',by.y='Sample.ID',all.x=TRUE)
dat$SUBTYPE = dat$BRCA_Subtype_PAM50
dat = dat[,-ncol(dat)]
clinSUB[['brca']] = dat


#### ACC ####
dat = clin[['acc']]
# subtypes
sub = SUBs[['acc']]
colnames(sub) = sub[6,]
sub = sub[-c(1:6),]
sub = subset(sub,select=c(SAMPLE,MethyLevel))
# merge
dat = merge(dat,sub,by.x='PATIENT_ID',by.y='SAMPLE',all.x=TRUE)
dat$SUBTYPE = dat$MethyLevel
dat = dat[,-ncol(dat)]
clinSUB[['acc']] = dat


#### BLCA ####
dat = clin[['blca']]
# subtypes
sub = SUBs[['blca']]
sub = subset(sub,select=c(Case.ID,mRNA.cluster))
# merge
dat = merge(dat,sub,by.x='PATIENT_ID',by.y='Case.ID',all.x=TRUE)
dat$SUBTYPE = dat$mRNA.cluster
dat = dat[,-ncol(dat)]
clinSUB[['blca']] = dat


#### CESC ####
dat = clin[['cesc']]
# subtype
sub = SUBs[['cesc']]
sub$X = gsub('-01A|-01B|-01C','',sub$X)
sub = subset(sub,select=c(X,SAMP.iCluster_All_k3))
# merge
dat = merge(dat,sub,by.x='PATIENT_ID',by.y='X',all.x=TRUE)
dat$SUBTYPE = ifelse(is.na(dat$SAMP.iCluster_All_k3),dat$SAMP.iCluster_All_k3,paste0('iCluster ',dat$SAMP.iCluster_All_k3))
dat = dat[,-ncol(dat)]
clinSUB[['cesc']] = dat


#### CHOL ####
dat = clin[['chol']]
# subtype
sub = SUBs[['chol']]
colnames(sub) = sub[1,]
sub = sub[-1,]
sub = subset(sub,select=c('Sample','mRNA Yulia -liver TumorMap k=7'))
# merge
dat = merge(dat,sub,by.x='PATIENT_ID',by.y='Sample',all.x=TRUE)
dat$SUBTYPE = ifelse(is.na(dat$`mRNA Yulia -liver TumorMap k=7`),dat$`mRNA Yulia -liver TumorMap k=7`,paste0('mRNA Cluster ',dat$`mRNA Yulia -liver TumorMap k=7`))
dat = dat[,-ncol(dat)]
clinSUB[['chol']] = dat


#### COADREAD ####
dat = clin[['coadread']]
# subtype
sub = SUBs[['esca_coad_stad_read']]
colnames(sub) = sub[1,]
sub = sub[-1,]
sub = sub[sub$`TCGA Project Code`=='COAD' | sub$`TCGA Project Code`=='READ',]
sub = subset(sub,select=c('TCGA Participant Barcode','Colorectal CMS'))
# merge
dat = merge(dat,sub,by.x='PATIENT_ID',by.y='TCGA Participant Barcode',all.x=TRUE)
dat$SUBTYPE = dat$`Colorectal CMS`
dat = dat[,-ncol(dat)]
clinSUB[['coadread']] = dat


#### DLBC ####  
clinSUB[['dlbc']] = clin[['dlbc']]


#### ESCA ####
dat = clin[['esca']]
# sub
sub = SUBs[['esca_coad_stad_read']]
colnames(sub) = sub[1,]
sub = sub[-1,]
sub = sub[sub$`TCGA Project Code`=='ESCA',]
sub = subset(sub,select=c('TCGA Participant Barcode','Hypermethylation category'))
# merge
dat = merge(dat,sub,by.x='PATIENT_ID',by.y='TCGA Participant Barcode',all.x=TRUE)
dat$SUBTYPE = dat$`Hypermethylation category`
dat = dat[,-ncol(dat)]
clinSUB[['esca']] = dat


#### HNSC ####
dat = clin[['hnsc']]
# subs
sub = SUBs[['hnsc']]
sub = subset(sub,select=c(Barcode,RNA))
sub$Barcode = gsub('\\.','-',sub$Barcode)
# merge
dat = merge(dat,sub,by.x='PATIENT_ID',by.y='Barcode',all.x=TRUE)
dat$SUBTYPE = dat$RNA
dat = dat[,-ncol(dat)]
clinSUB[['hnsc']] = dat


#### KICH ####
dat = clin[['kich']]
# subs
sub = SUBs[['kirc_kirp_kich']]
sub = sub[sub$Original.TCGA.project=='KICH',]
sub = subset(sub,select=c(bcr_patient_barcode,mRNA.Cluster))
# merge
dat = merge(dat,sub,by.x='PATIENT_ID',by.y='bcr_patient_barcode',all.x=TRUE)
dat$SUBTYPE = ifelse(is.na(dat$mRNA.Cluster),dat$mRNA.Cluster,paste0('mRNA Cluster ',dat$mRNA.Cluster)) 
dat = dat[,-ncol(dat)]
clinSUB[['kich']] = dat


#### KIRC ####
dat = clin[['kirc']]
# subtypes
sub = SUBs[['kirc_kirp_kich']]
sub = sub[sub$Original.TCGA.project=='KIRC',]
sub = subset(sub,select=c(bcr_patient_barcode,mRNA.Cluster))
# merge
dat = merge(dat,sub,by.x='PATIENT_ID',by.y='bcr_patient_barcode',all.x=TRUE)
dat$mRNA.Cluster[dat$mRNA.Cluster=='n/a']=NA
dat$SUBTYPE = ifelse(is.na(dat$mRNA.Cluster),dat$mRNA.Cluster,paste0('mRNA Cluster ',dat$mRNA.Cluster)) 
dat = dat[,-ncol(dat)]
clinSUB[['kirc']] = dat


#### KIRP ####
dat = clin[['kirp']]
# subtypes
sub = SUBs[['kirc_kirp_kich']]
sub = sub[sub$Original.TCGA.project=='KIRP',]
sub = subset(sub,select=c(bcr_patient_barcode,mRNA.Cluster))
# merge
dat = merge(dat,sub,by.x='PATIENT_ID',by.y='bcr_patient_barcode',all.x=TRUE)
dat$mRNA.Cluster[dat$mRNA.Cluster=='n/a']=NA
dat$SUBTYPE = ifelse(is.na(dat$mRNA.Cluster),dat$mRNA.Cluster,paste0('mRNA Cluster ',dat$mRNA.Cluster)) 
dat = dat[,-ncol(dat)]
clinSUB[['kirp']] = dat


#### LAML ####
dat = clin[['laml']]
# sub
sub = SUBs[['laml']]
sub$TCGA.Patient.ID = paste0('TCGA-AB-',sub$TCGA.Patient.ID)
sub = subset(sub,select=c(TCGA.Patient.ID,Molecular.Classification))
# merge
dat = merge(dat,sub,by.x='PATIENT_ID',by.y='TCGA.Patient.ID',all.x=TRUE)
dat$SUBTYPE = dat$Molecular.Classification
dat = dat[,-ncol(dat)]
clinSUB[['laml']] = dat


#### LGG ####
dat = clin[['lgg']]
# sub
sub = SUBs[['gbm_lgg']]
colnames(sub) = sub[1,]
sub = sub[-1,]
sub = sub[sub$Study=='Brain Lower Grade Glioma',]
sub = subset(sub,select=c('Case','Transcriptome Subtype'))
# merge
dat = merge(dat,sub,by.x='PATIENT_ID',by.y='Case',all.x=TRUE)
dat$SUBTYPE = dat$`Transcriptome Subtype`
dat = dat[,-ncol(dat)]
clinSUB[['lgg']] = dat


#### LIHC ####
dat = clin[['lihc']]
# sub
sub = SUBs[['lihc']]
colnames(sub) = sub[3,]
sub = sub[-c(1:3),]
sub$Barcode = gsub('-01A|-01B','',sub$Barcode)
sub = subset(sub,select=c('Barcode','iCluster clusters (k=3, Ronglai Shen)'))
# merge
dat = merge(dat,sub,by.x='PATIENT_ID','Barcode',all.x=TRUE)
dat$SUBTYPE = dat$`iCluster clusters (k=3, Ronglai Shen)`
dat = dat[,-ncol(dat)]
clinSUB[['lihc']] = dat


#### LUAD ####
dat = clin[['luad']]
# sub
sub = SUBs[['luad']]
colnames(sub) = sub[4,]
sub = sub[-c(1:4),]
sub = subset(sub,select=c('Tumor ID','expression_subtype'))
# merge
dat = merge(dat,sub,by.x='PATIENT_ID',by.y='Tumor ID',all.x=TRUE)
dat$SUBTYPE = dat$expression_subtype
dat = dat[,-ncol(dat)]
clinSUB[['luad']] = dat


#### LUSC ####
dat = clin[['lusc']]
# sub
sub = SUBs[['lusc']]
colnames(sub) = sub[2,]
sub = sub[-c(1:2),]
sub$`Tumor ID` = gsub('LUSC','TCGA',sub$`Tumor ID`)
colnames(sub)[101] = 'sub'
sub = subset(sub,select=c('Tumor ID','sub'))
# merge
dat = merge(dat,sub,by.x='PATIENT_ID',by.y='Tumor ID',all.x=TRUE)
dat$SUBTYPE = dat$sub
dat = dat[,-ncol(dat)]
clinSUB[['lusc']] = dat


#### MESO ####
dat = clin[['meso']]
# sub
sub = SUBs[['meso']]
sub$TCGA_barcode = gsub('-01','',sub$TCGA_barcode)
sub = subset(sub,select=c(TCGA_barcode,iCluster_k4.5types))
# merge
dat = merge(dat,sub,by.x='PATIENT_ID',by.y='TCGA_barcode',all.x=TRUE)
dat$SUBTYPE = ifelse(is.na(dat$iCluster_k4.5types),dat$iCluster_k4.5types,paste0('iCluster ',dat$iCluster_k4.5types))  
dat = dat[,-ncol(dat)]
clinSUB[['meso']] = dat


#### OV ####
dat = clin[['ov']]
# sub
sub = SUBs[['brca_cesc_ov_ucec_ucs']]
colnames(sub) = sub[1,]
sub = sub[-1,]
sub = subset(sub,select=c(Sample.ID,OV_Subtype))
# merge
dat = merge(dat,sub,by.x='PATIENT_ID',by.y='Sample.ID',all.x=TRUE)
dat$SUBTYPE = dat$OV_Subtype
dat = dat[,-ncol(dat)]
clinSUB[['ov']] = dat


#### PAAD ####
dat = clin[['paad']]
# sub
sub = SUBs[['paad']]
colnames(sub) = sub[1,]
sub = sub[-1,]
sub$`Tumor Sample ID` = gsub('-01A','',sub$`Tumor Sample ID`)
sub = subset(sub,select=c('Tumor Sample ID','mRNA Bailey Clusters (All 150 Samples) 1squamous 2immunogenic 3progenitor 4ADEX'))
# merge
dat = merge(dat,sub,by.x='PATIENT_ID',by.y='Tumor Sample ID',all.x=TRUE)
dat$SUBTYPE = ifelse(is.na(dat$`mRNA Bailey Clusters (All 150 Samples) 1squamous 2immunogenic 3progenitor 4ADEX`),dat$`mRNA Bailey Clusters (All 150 Samples) 1squamous 2immunogenic 3progenitor 4ADEX`,paste0('mRNA Bailey Cluster ',dat$`mRNA Bailey Clusters (All 150 Samples) 1squamous 2immunogenic 3progenitor 4ADEX`)) 
dat = dat[,-ncol(dat)]
clinSUB[['paad']] = dat


#### PCPG ####
dat = clin[['pcpg']]
# sub
sub = SUBs[['pcpg']]
colnames(sub) = sub[2,]
sub = sub[-c(1:2),]
sub$`Sample ID` = gsub('-01A|-01B|-05A|-06A','',sub$`Sample ID`)
sub = subset(sub,select=c('Sample ID','mRNA Subtype Clusters'))
# merge
dat = merge(dat,sub,by.x='PATIENT_ID',by.y='Sample ID',all.x=TRUE)
dat$SUBTYPE = dat$`mRNA Subtype Clusters`
dat = dat[,-ncol(dat)]
clinSUB[['pcpg']] = dat


#### PRAD ####
dat = clin[['prad']]
#
sub = SUBs[['prad']]
sub = subset(sub,select=c(PATIENT_ID,Subtype))
# merge
dat = merge(dat,sub,by='PATIENT_ID',all.x=TRUE)
dat$SUBTYPE = dat$Subtype
dat = dat[,-ncol(dat)]
clinSUB[['prad']] = dat


#### SARC ####
dat = clin[['sarc']]
# 
sub = SUBs[['sarc']]
colnames(sub) = sub[1,]
sub = sub[-1,]
sub$`TCGA barcode` = gsub('-01','',sub$`TCGA barcode`)
sub = subset(sub,select=c('TCGA barcode','mRNA cluster'))
# merge
dat = merge(dat,sub,by.x='PATIENT_ID',by.y='TCGA barcode',all.x=TRUE)
dat$SUBTYPE = ifelse(is.na(dat$`mRNA cluster`),dat$`mRNA cluster`,paste0('mRNA Cluster ',dat$`mRNA cluster`)) 
dat = dat[,-ncol(dat)]
clinSUB[['sarc']] = dat


#### SKCM ####
dat = clin[['skcm']]
# sub
sub = SUBs[['skcm']]
colnames(sub) = sub[1,]
sub = sub[-1,]
sub$Name = gsub('-01|-06','',sub$Name)
sub = subset(sub,select=c(Name,MUTATIONSUBTYPES))
# merge
dat = merge(dat,sub,by.x='PATIENT_ID',by.y='Name',all.x=TRUE)
dat$SUBTYPE = dat$MUTATIONSUBTYPES
dat = dat[,-ncol(dat)]
clinSUB[['skcm']] = dat


#### STAD ####
dat = clin[['stad']]
# subs
sub = SUBs[['esca_coad_stad_read']]
colnames(sub) = sub[1,]
sub = sub[-1,]
sub = sub[sub$`TCGA Project Code`=='STAD',]
sub = subset(sub,select=c('TCGA Participant Barcode','Molecular_Subtype'))
# merge
dat = merge(dat,sub,by.x='PATIENT_ID',by.y='TCGA Participant Barcode',all.x=TRUE)
dat$SUBTYPE = dat$Molecular_Subtype
dat = dat[,-ncol(dat)]
clinSUB[['stad']] = dat


#### TGCT ####
dat = clin[['tgct']]
# sub
sub = SUBs[['tgct']]
sub$sample = gsub('-01|-05','',sub$sample)
sub = subset(sub,select=c(sample,methylation_k5))
# merge
dat = merge(dat,sub,by.x='PATIENT_ID',by.y='sample',all.x = TRUE)
dat$SUBTYPE = dat$methylation_k5
dat = dat[,-ncol(dat)]
clinSUB[['tgct']] = dat


#### THCA ####
dat = clin[['thca']]
# sub
sub = SUBs[['thca']]
sub = subset(sub,select=c(sample,mRNA_Cluster_number))
sub$sample = gsub('-01A|-01B','',sub$sample)
# merge
dat = merge(dat,sub,by.x='PATIENT_ID',by.y='sample',all.x=TRUE)
dat$SUBTYPE = ifelse(is.na(dat$mRNA_Cluster_number),dat$mRNA_Cluster_number,paste0('mRNA Cluster ',dat$mRNA_Cluster_number)) 
dat = dat[,-ncol(dat)]
clinSUB[['thca']] = dat


#### THYM ####
dat = clin[['thym']]
# sub
sub = SUBs[['thym']]
colnames(sub) = sub[1,]
sub = sub[-1,]
sub$`Tumor Barcode` = gsub('-01A','',sub$`Tumor Barcode`)
sub = subset(sub,select=c('Tumor Barcode','mRNA clusters (4)'))
# merge
dat = merge(dat,sub,by.x='PATIENT_ID',by.y='Tumor Barcode',all.x=TRUE)
dat$SUBTYPE = ifelse(is.na(dat$`mRNA clusters (4)`),dat$`mRNA clusters (4)`,paste0('mRNA Cluster ',dat$`mRNA clusters (4)`)) 
dat = dat[,-ncol(dat)]
clinSUB[['thym']] = dat


#### UCEC ####
dat = clin[['ucec']]
# sub
sub = SUBs[['ucec']]
sub = subset(sub,select=c(bcr_patient_barcode,IntegrativeCluster))
# merge
dat = merge(dat,sub,by.x='PATIENT_ID','bcr_patient_barcode',all.x=TRUE)
dat$SUBTYPE = dat$IntegrativeCluster
dat = dat[,-ncol(dat)]
clinSUB[['ucec']] = dat


#### UCS ####
dat = clin[['ucs']]
# 
sub = SUBs[['ucs']]
colnames(sub) = sub[1,]
sub = sub[-1,]
sub = subset(sub,select=c('bcr_patient_barcode','Histologic classification'))
# merge
dat = merge(dat,sub,by.x='PATIENT_ID','bcr_patient_barcode',all.x=TRUE)
dat$SUBTYPE = dat$`Histologic classification`
dat = dat[,-ncol(dat)]
clinSUB[['ucs']] = dat


#### UVM ####
dat = clin[['uvm']]
# sub
sub = SUBs[['uvm']]
colnames(sub) = sub[1,]
sub = sub[-1,]
sub = subset(sub,select=c('Patient ID','SCNA Cluster No.'))
# merge
dat = merge(dat,sub,by.x='PATIENT_ID',by.y='Patient ID',all.x=TRUE)
dat$SUBTYPE = ifelse(is.na(dat$`SCNA Cluster No.`),dat$`SCNA Cluster No.`,paste0('SCNA Cluster ',dat$`SCNA Cluster No.`)) 
dat = dat[,-ncol(dat)]
clinSUB[['uvm']] = dat

# put dfs in alphabetical order
clinSUB = clinSUB[order(names(clinSUB))]



#### Add metastatic / primary info from SAMPLE SPECIFIC csvs ####
# read in SAMPLE SPECIFIC clinical data
path = '/Volumes/BioDiscovery/Kristin/SurVarCN_previousVersions/SurVarCN_updatedTCGA/cBioPortal'
temp = list.files(path = path,pattern = '*data_bcr_clinical_data_sample.txt',full.names = TRUE)
SAMPLES = lapply(temp, read.delim, stringsAsFactors = F, sep='\t')
names(SAMPLES) = gsub('_tcga_data_bcr_clinical_data_sample.txt','',basename(temp))
#restructure df
for (i in 1:length(SAMPLES)) {
  colnames(SAMPLES[[i]]) = SAMPLES[[i]][4,]
  SAMPLES[[i]] = SAMPLES[[i]][-c(1:4),]
  SAMPLES[[i]] = subset(SAMPLES[[i]],select=c(PATIENT_ID,SAMPLE_ID,SAMPLE_TYPE))
}

# table number of patients with metastatic vs primary samples
# tab_meta = data.frame(matrix(ncol=3,nrow=length(SAMPLES)))
# colnames(tab_meta) = c('Cancer','Metastasis','Primary')
# for (i in 1:length(SAMPLES)) {
#   tab_meta[i,1] = toupper(names(SAMPLES)[i])
#   tab_meta[i,2] = length(SAMPLES[[i]]$SAMPLE_TYPE[SAMPLES[[i]]$SAMPLE_TYPE=='Metastasis'])
#   tab_meta[i,3] = length(SAMPLES[[i]]$SAMPLE_TYPE[SAMPLES[[i]]$SAMPLE_TYPE=='Primary'])
# }
# write.csv(tab_meta,'/Volumes/BioDiscovery/Kristin/SurVarCN_Metastatic/metaPrim_table.csv',row.names = F)


# Will only work with heatmaps, too little metastatic samples (except skcm) for survival plots
for (i in 1:length(clinSUB)) {
  clinSUB[[i]] = merge(clinSUB[[i]],SAMPLES[[i]],by='PATIENT_ID',all.x=TRUE)
  clinSUB[[i]] = clinSUB[[i]][,c(15,2:10,16,11:14)]
}



#### Get clinical RT info from RTCGA package #####
library(RTCGA)
library(RTCGA.clinical)

# download clinical data from package environment
clinCo = list()
for (obj in ls(package:RTCGA.clinical)) {
  clinCo[[obj]] = get(obj)
}
# change name to TCGA cancer id
# subset patient id and RT columns
# if RT column non existent, create one and fill with NA
for (i in 1:length(clinCo)) {
  names(clinCo)[i] = tolower(gsub('.clinical','',names(clinCo)[i]))
  if ('patient.follow_ups.follow_up.radiation_therapy' %in% colnames(clinCo[[i]])) {
    clinCo[[i]] = subset(clinCo[[i]],select=c(patient.bcr_patient_barcode,patient.follow_ups.follow_up.radiation_therapy))
  } else {
    clinCo[[i]]$patient.follow_ups.follow_up.radiation_therapy = NA
    clinCo[[i]] = subset(clinCo[[i]],select=c(patient.bcr_patient_barcode,patient.follow_ups.follow_up.radiation_therapy))
  }
  clinCo[[i]]$patient.bcr_patient_barcode = toupper(clinCo[[i]]$patient.bcr_patient_barcode)
}
# subset only cancers in our dataset
clinCo = clinCo[names(clinCo) %in% names(clinSUB)]
#
# Check #NAs in original and new data, replace RT if nec
for (i in 1:length(clinSUB)) {
  origNA = length(clinSUB[[i]]$RADIATION_TREATMENT_ADJUVANT[is.na(clinSUB[[i]]$RADIATION_TREATMENT_ADJUVANT)])
  newNA = length(clinCo[[i]]$patient.follow_ups.follow_up.radiation_therapy[is.na(clinCo[[i]]$patient.follow_ups.follow_up.radiation_therapy)])
  if(origNA > newNA) {
    clinSUB[[i]]$temp_sample_id = gsub('.{3}$', '', clinSUB[[i]]$SAMPLE_ID)
    clinSUB[[i]] = merge(clinSUB[[i]],clinCo[[i]],by.x='temp_sample_id',by.y='patient.bcr_patient_barcode',all.x = TRUE)
    clinSUB[[i]]$RADIATION_TREATMENT_ADJUVANT = clinSUB[[i]]$patient.follow_ups.follow_up.radiation_therapy
    clinSUB[[i]] = subset(clinSUB[[i]],select=-c(patient.follow_ups.follow_up.radiation_therapy,temp_sample_id))
    clinSUB[[i]]$RADIATION_TREATMENT_ADJUVANT = toupper(clinSUB[[i]]$RADIATION_TREATMENT_ADJUVANT)
  }
}


#### Plot clinical info for each cancer ####

#### Histograms ####

relevantClin = clinSUB
for (i in 1:length(clinSUB)) {
  relevantClin[[i]] = relevantClin[[i]][,c(2:12,14)]
  notEmpty = apply(relevantClin[[i]], 2, function(x) !all(is.na(x)))
  relevantClin[[i]] = relevantClin[[i]][,notEmpty]
  
  for (k in 1:ncol(relevantClin[[i]])) {
    relevantClin[[i]][,k] =  str_trunc(relevantClin[[i]][,k], 15, "center",ellipsis = '..')
  }
}
tableClinical = list()
pl = list()
grids = list()
grids_display = list()
colors = primary.colors(12)
pdf('/Volumes/BioDiscovery/Kristin/SurVarCN_ClinicalData/TCGA_Clinical_Histograms.pdf')
for (i in 1:length(relevantClin)) {
  tableClinical[[i]] = list()
  pl[[i]] = list()
  for (j in 1:length(colnames(relevantClin[[i]]))) {
    tableClinical[[i]][[j]] = data.frame(table(relevantClin[[i]][j],useNA='ifany'))
    names(tableClinical[[i]])[j] = colnames(relevantClin[[i]])[j]
    
    pl[[i]][[j]] = ggplot(tableClinical[[i]][[j]]) + 
      geom_bar(aes(x=Var1,y=Freq), stat = 'identity',fill = colors[j]) + 
      theme_bw() +
      xlab(names(tableClinical[[i]][j])) +
      ylab('# Patients') +
      theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, size = 6, hjust = 1), axis.title.x = element_text(size = 8))
    names(pl[[i]])[j] = colnames(relevantClin[[i]])[j]
  }
  # grids[[i]] = arrangeGrob(grobs=pl[[i]],top=toupper(names(clinSUB)[i]))
  # grids_display[[i]] = grid.arrange(grobs=pl[[i]],top=toupper(names(clinSUB)[i]))
  # ggsave(file=paste0('/Volumes/BioDiscovery/Kristin/SurVarCN_ClinicalData/',toupper(names(clinSUB)[i]),'.pdf'),grids[[i]])
  
  grid.arrange(grobs=pl[[i]],top=toupper(names(clinSUB)[i]))
}
dev.off()
names(tableClinical) = toupper(names(clinSUB))
names(pl) = toupper(names(clinSUB))
# names(grids) = toupper(names(clinSUB))



addSmallLegend <- function(myPlot, pointSize = 2, textSize = 6, spaceLegend = 0.5) {
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_text(size = textSize), 
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"))
}




#### Pie charts ####
relevantClin = clinSUB
for (i in 1:length(clinSUB)) {
  relevantClin[[i]] = relevantClin[[i]][,c(2:12,14)]
  notEmpty = apply(relevantClin[[i]], 2, function(x) !all(is.na(x)))
  relevantClin[[i]] = relevantClin[[i]][,notEmpty]
  
  for (k in 1:ncol(relevantClin[[i]])) {
    relevantClin[[i]][,k] =  str_trunc(relevantClin[[i]][,k], 30, "center",ellipsis = '..')
  }
}
tableClinical = list()
pl = list()
grids = list()
grids_display = list()
colors = primary.colors(12)
pdf('/Volumes/BioDiscovery/Kristin/SurVarCN_ClinicalData/TCGA_Clinical_PieCharts.pdf')
for (i in 1:length(relevantClin)) {
  tableClinical[[i]] = list()
  pl[[i]] = list()
  for (j in 1:length(colnames(relevantClin[[i]]))) {
    tableClinical[[i]][[j]] = data.frame(table(relevantClin[[i]][j],useNA='ifany'))
    names(tableClinical[[i]])[j] = colnames(relevantClin[[i]])[j]
    # colnames(tableClinical[[i]][[j]])[1] = colnames(relevantClin[[i]])[j]
    
    pl[[i]][[j]] = ggplot(tableClinical[[i]][[j]],aes(x='',y=Freq,fill=Var1)) + 
      geom_bar(width = 1 ,stat = 'identity') + 
      theme_bw() +
      theme(axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text = element_blank()) +
      coord_polar('y') +
      guides(fill=guide_legend(title=colnames(relevantClin[[i]])[j]))
      # geom_text(aes(label = paste(round(Freq / sum(Freq) * 100, 1), "%")),
      #           position = position_stack(vjust = 0.5), size = 2)
    pl[[i]][[j]] = addSmallLegend(pl[[i]][[j]])
    names(pl[[i]])[j] = colnames(relevantClin[[i]])[j]
  }
  p = plot_grid(plotlist = pl[[i]],align = 'v',ncol = 2)
  # now add the title
  title <- ggdraw() + 
    draw_label(toupper(names(clinSUB)[i])
      ,
      fontface = 'bold',
      x = 0,
      hjust = 0
    ) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 7)
    )
  
  print(plot_grid(
    title, p,
    ncol = 1,
    # rel_heights values control vertical title margins
    rel_heights = c(0.1, 1)
  ))
  # grid.arrange(grobs=pl[[i]],top=toupper(names(clinSUB)[i]),ncol=2)
}
dev.off()
names(tableClinical) = toupper(names(clinSUB))
names(pl) = toupper(names(clinSUB))




