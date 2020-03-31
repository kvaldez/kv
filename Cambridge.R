library(Hmisc)
library(rlist)
library(tidyverse)
library(purrr)
library(amap)
library(pheatmap)
library(plotly)
library(ggplot2)
library(reshape2)
library(dplyr)

#### Read in mafs - Vaf heatmap - only for ALL non-filtered variants ####
temp = list.files(path = '/Users/valdezkm/Documents/Cambridge/mutect2',pattern = '*.maf',full.names = TRUE)
allMafs = lapply(temp, read.delim, stringsAsFactors = F)
names(allMafs) = basename(temp)
pheno = data.frame(fileNames = basename(temp))
pheno$fileNames = gsub('.*\\+|_.*','',pheno$fileNames)
seqSampleKey = read.delim('/Users/valdezkm/Documents/Cambridge/SequencingSampleKey_McAbeeJH.txt')
pheno = merge(pheno,seqSampleKey,by.x='fileNames',by.y='File.Letter.Name',all.x=TRUE, sort=FALSE)
pheno$Label.Name = gsub('DNA|\\*|12/16|12/2016','',pheno$Label.Name)
pheno$Label.Name = gsub(' ','',pheno$Label.Name)
pheno$Label.Name = gsub('_\\(2\\)','',pheno$Label.Name)

# make sure all above worked right:
# pheno$fileNames == gsub('.*\\+|_.*','',names(allMafs))

# rename maf names
names(allMafs) = pheno$Label.Name

# GBM driver genes
# JoeGenes = read.delim('/Users/valdezkm/Documents/Driver_Genes/Joseph_GBM_Drivers.txt')
# JoeGenes = as.character(unique(JoeGenes$GBM.Drive.Genes))



#### Write out UNFILTERED variants (no germline) ####
for (i in 1:length(allMafs)) {
  write.csv(allMafs[[i]],paste0('/Volumes/BioDiscovery/Lab/Tofilon/P06_Cambridge_Joseph/mutect2/UNfiltered_variants/',names(allMafs)[i],'_maf.csv'),row.names = F)
}



# add vaf
for (i in 1:length(allMafs)) {
  allMafs[[i]]$vaf = allMafs[[i]]$t_alt_count / allMafs[[i]]$t_depth
  allMafs[[i]] = allMafs[[i]][allMafs[[i]]$vaf > 0.05,]
}

#### Pearsons correlation coefficients - VAFs - level 1 filtered data (germline variants removed) ####
justVAFs = c()
for (i in 1:length(allMafs)) {
  justVAFs[[i]] = subset(allMafs[[i]],select=c(Hugo_Symbol,Chromosome,Start_Position,vaf))
  colnames(justVAFs[[i]])[4] = names(allMafs[i])
}
df_justVAFs = justVAFs %>% reduce(full_join, by = c('Hugo_Symbol','Chromosome','Start_Position'))
mat_justVAFs = as.matrix(df_justVAFs[,4:ncol(df_justVAFs)])

# jaccard distance
jac_mat = mat_justVAFs
jac_mat[is.na(jac_mat)] = 0         # need to do this for jaccard binary correlation
jac = dist(t(jac_mat),method='binary',diag = TRUE,upper = TRUE)
jac = as.matrix(jac)
jac = 1-jac

pheatmap(jac,display_numbers = T,fontsize = 14)


# J3 VAF heatmap - pearson
J3 = justVAFs[grepl('J3',names(allMafs))]
df_justVAFs = J3 %>% reduce(full_join, by = c('Hugo_Symbol','Chromosome','Start_Position'))
mat_justVAFs = as.matrix(df_justVAFs[,4:ncol(df_justVAFs)])
cors = rcorr(mat_justVAFs, type = 'pearson')
cors_r = cors$r
pheatmap(cors_r,display_numbers = T,fontsize = 14)
#
# J3 jaccard
# mat_justVAFs[is.na(mat_justVAFs)] = 0
# jac = dist(t(mat_justVAFs),method='binary',diag=TRUE,upper=TRUE)
# jac = as.matrix(jac)
# jac = 1-jac
# pheatmap(jac,display_numbers = T,fontsize = 14)

#
# J7 VAF heatmap
J7 = justVAFs[grepl('J7',names(allMafs))]
df_justVAFs = J7 %>% reduce(full_join, by = c('Hugo_Symbol','Chromosome','Start_Position'))
mat_justVAFs = as.matrix(df_justVAFs[,4:ncol(df_justVAFs)])
cors = rcorr(mat_justVAFs, type = 'pearson')
cors_r = cors$r
pheatmap(cors_r,display_numbers = T,fontsize = 14)
#
# J7 jaccard
mat_justVAFs[is.na(mat_justVAFs)] = 0
jac = dist(t(mat_justVAFs),method='binary',diag=TRUE,upper=TRUE)
jac = as.matrix(jac)
jac = 1-jac
pheatmap(jac,display_numbers = T,fontsize = 14)
#
# J14 VAF heatmap
J14 = justVAFs[grepl('J14',names(allMafs))]
df_justVAFs = J14 %>% reduce(full_join, by = c('Hugo_Symbol','Chromosome','Start_Position'))
mat_justVAFs = as.matrix(df_justVAFs[,4:ncol(df_justVAFs)])
cors = rcorr(mat_justVAFs, type = 'pearson')
cors_r = cors$r
pheatmap(cors_r,display_numbers = T,fontsize = 14)
#
# J14 jaccard
mat_justVAFs[is.na(mat_justVAFs)] = 0
jac = dist(t(mat_justVAFs),method='binary',diag=TRUE,upper=TRUE)
jac = as.matrix(jac)
jac = 1-jac
pheatmap(jac,display_numbers = T,fontsize = 14)




#### Pearsons correlation coefficients - VAFs - level 2 filtered data (germline and non-protein-altering variants removed) ####
justVAFs = c()
for (i in 1:length(filtMafs)) {
  justVAFs[[i]] = subset(filtMafs[[i]],select=c(Hugo_Symbol,Chromosome,Start_Position,vaf))
  colnames(justVAFs[[i]])[4] = names(filtMafs[i])
}
df_justVAFs = justVAFs %>% reduce(full_join, by = c('Hugo_Symbol','Chromosome','Start_Position'))
mat_justVAFs = as.matrix(df_justVAFs[,4:ncol(df_justVAFs)])

# jaccard distance
jac_mat = mat_justVAFs
jac_mat[is.na(jac_mat)] = 0         # need to do this for jaccard binary correlation
jac = dist(t(jac_mat),method='binary',diag = TRUE,upper = TRUE)
jac = as.matrix(jac)
jac = 1-jac
pheatmap(jac,display_numbers = T,fontsize = 14)


# J3 VAF heatmap
J3 = justVAFs[grepl('J3',names(filtMafs))]
df_justVAFs = J3 %>% reduce(full_join, by = c('Hugo_Symbol','Chromosome','Start_Position'))
mat_justVAFs = as.matrix(df_justVAFs[,4:ncol(df_justVAFs)])
cors = rcorr(mat_justVAFs, type = 'pearson')
cors_r = cors$r
pheatmap(cors_r,display_numbers = T,fontsize = 14)
#
# J3 jaccard
# mat_justVAFs[is.na(mat_justVAFs)] = 0
# jac = dist(t(mat_justVAFs),method='binary',diag=TRUE,upper=TRUE)
# jac = as.matrix(jac)
# jac = 1-jac
# pheatmap(jac,display_numbers = T,fontsize = 14)
#
# J7 VAF heatmap
J7 = justVAFs[grepl('J7',names(filtMafs))]
df_justVAFs = J7 %>% reduce(full_join, by = c('Hugo_Symbol','Chromosome','Start_Position'))
mat_justVAFs = as.matrix(df_justVAFs[,4:ncol(df_justVAFs)])
cors = rcorr(mat_justVAFs, type = 'pearson')
cors_r = cors$r
pheatmap(cors_r,display_numbers = T,fontsize = 14)
#
# J7 jaccard
# mat_justVAFs[is.na(mat_justVAFs)] = 0
# jac = dist(t(mat_justVAFs),method='binary',diag=TRUE,upper=TRUE)
# jac = as.matrix(jac)
# jac = 1-jac
# pheatmap(jac,display_numbers = T,fontsize = 14)
#
# J14 VAF heatmap
J14 = justVAFs[grepl('J14',names(filtMafs))]
df_justVAFs = J14 %>% reduce(full_join, by = c('Hugo_Symbol','Chromosome','Start_Position'))
mat_justVAFs = as.matrix(df_justVAFs[,4:ncol(df_justVAFs)])
cors = rcorr(mat_justVAFs, type = 'pearson')
cors_r = cors$r
pheatmap(cors_r,display_numbers = T,fontsize = 14)
#
# J14 jaccard
# mat_justVAFs[is.na(mat_justVAFs)] = 0
# jac = dist(t(mat_justVAFs),method='binary',diag=TRUE,upper=TRUE)
# jac = as.matrix(jac)
# jac = 1-jac
# pheatmap(jac,display_numbers = T,fontsize = 14)







#### 3D PCA ####
tedf= t(mat_justVAFs)

# removes zero  variances (issue with small sample sizes)
if (length(which(apply(tedf, 2, var)==0)) >= 0){
  tedf = tedf[ , apply(tedf, 2, var) != 0]
}

pca=prcomp(tedf, scale. = T)

pc1.var=round(pca$sdev[1]^2/sum(pca$sdev^2)*100,2)
pc2.var=round(pca$sdev[2]^2/sum(pca$sdev^2)*100,2)
pc3.var=round(pca$sdev[3]^2/sum(pca$sdev^2)*100,2)

group.v=as.vector(colnames(mat_justVAFs))
p = plotly::plot_ly(as.data.frame(pca$x[,1:3]), x = ~PC1, y = ~PC2, z = ~PC3, hoverinfo="text",
                    hovertext = ~group.v) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = paste0("PC1 (",pc1.var,"%)")),
                      yaxis = list(title = paste0("PC2 (",pc2.var,"%)")),
                      zaxis = list(title = paste0("PC3 (",pc3.var,"%)"))))
p


#### Filtered MAFs -  script is same as allMafs but includes filtering -> Heatmap of GBM driver genes -####
#### Start here for YAPSA prep ####
temp = list.files(path = '/Users/valdezkm/Documents/Cambridge/mutect2',pattern = '*.maf',full.names = TRUE)
gbmMafs = lapply(temp, read.delim, stringsAsFactors = F)
names(gbmMafs) = basename(temp)
pheno = data.frame(fileNames = basename(temp))
pheno$fileNames = gsub('.*\\+|_.*','',pheno$fileNames)
seqSampleKey = read.delim('/Users/valdezkm/Documents/Cambridge/SequencingSampleKey_McAbeeJH.txt')
pheno = merge(pheno,seqSampleKey,by.x='fileNames',by.y='File.Letter.Name',all.x=TRUE, sort=FALSE)
pheno$Label.Name = gsub('DNA|\\*|12/16|12/2016','',pheno$Label.Name)
pheno$Label.Name = gsub(' ','',pheno$Label.Name)
pheno$Label.Name = gsub('_\\(2\\)','',pheno$Label.Name)
#
# make sure all above worked right:
# pheno$fileNames == gsub('.*\\+|_.*','',names(gbmMafs))
#
# rename maf names
names(gbmMafs) = pheno$Label.Name
#
# mean of num mutations in each sample
num = lapply(gbmMafs, function(x) nrow(x))
mean(unlist(num))
#
# add vaf, filter out low vafs
for (i in 1:length(gbmMafs)) {
  gbmMafs[[i]]$vaf = gbmMafs[[i]]$t_alt_count / gbmMafs[[i]]$t_depth
  gbmMafs[[i]] = gbmMafs[[i]][gbmMafs[[i]]$vaf > 0.05,]
}
#
# filter out non-protein-altering variants
variantsToKeep = c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site","Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del","In_Frame_Ins", "Missense_Mutation")
filtMafs = list()
for (i in 1:length(gbmMafs)) {
  filtMafs[[i]] = gbmMafs[[i]][gbmMafs[[i]]$Variant_Classification %in% variantsToKeep,]
  names(filtMafs)[i] = names(gbmMafs)[i]
}
#### End filtering mafs ####
#### END YAPSA PREP ####


#### Write out filtered variants (no germline or non-protein-altering) ####
colsToKeep = c(1:19,35:76,133)
temp = filtMafs
for (i in 1:length(filtMafs)) {
  temp[[i]] = temp[[i]][,colsToKeep]
  write.csv(temp[[i]],paste0('/Users/valdezkm/Documents/Cambridge/mutect2/filtered_variants/',names(temp)[i],'_filteredVariants.csv'),row.names = F)
}


# mean of num mutations in each sample
num = lapply(filtMafs, function(x) nrow(x))
mean(unlist(num))

# GBM driver genes
# JoeGenes = read.delim('/Users/valdezkm/Documents/Driver_Genes/Joseph_GBM_Drivers.txt')
# JoeGenes = as.character(unique(JoeGenes$GBM.Drive.Genes))
# # subset driver genes
# for (i in 1:length(filtMafs)) {
#   filtMafs[[i]] = filtMafs[[i]][filtMafs[[i]]$Hugo_Symbol %in% JoeGenes,]
# }

# heatmap of GBM variants and vafs
GBMvar = list()
for (i in 1:length(filtMafs)) {
  if (nrow(filtMafs[[i]]) > 0) {
    GBMvar[[i]] = subset(filtMafs[[i]],select=c(Hugo_Symbol,HGVSc,vaf))
    GBMvar[[i]]$geneVar = paste(GBMvar[[i]]$Hugo_Symbol,GBMvar[[i]]$HGVSc)
    GBMvar[[i]]$Sample = names(filtMafs)[i]
  } else {
    GBMvar[[i]] = data.frame(matrix(ncol=3))
    colnames(GBMvar[[i]]) = c('geneVar','vaf')
    GBMvar[[i]]$geneVar = as.character(GBMvar[[i]]$geneVar)
    GBMvar[[i]]$vaf = as.numeric(GBMvar[[i]]$vaf)
    GBMvar[[i]]$Sample = names(filtMafs)[i]
  }
  GBMvar[[i]] = subset(GBMvar[[i]],select=c(geneVar,vaf,Sample))
  names(GBMvar)[i] = names(filtMafs)[i]
}

df_GBMvar = GBMvar %>% reduce(full_join,by=c('geneVar','vaf','Sample'))
df_GBMvar = dcast(df_GBMvar,geneVar~Sample,value.var = 'vaf')
df_GBMvar = df_GBMvar[!is.na(df_GBMvar$geneVar),]
rownames(df_GBMvar) = df_GBMvar$geneVar
df_GBMvar = df_GBMvar[,-1]
textVar = df_GBMvar
textVar = round(textVar,digits = 2)
textVar[is.na(textVar)] = ''

pheatmap(df_GBMvar,treeheight_row = 0,treeheight_col = 0,display_numbers = textVar, fontsize = 12, cluster_cols = F,cluster_rows = FALSE, na_col = 'white')


# tables of variants per sample (unfiltered)
numRows = lapply(gbmMafs,function(x) {unlist(table(x$Variant_Classification))})
numRows = max(unlist(lapply(numRows, nrow)))
allVars_table = data.frame(matrix(ncol=12,nrow=numRows))

allVars_table = list()
for (i in 1:length(gbmMafs)) {
  allVars_table[[i]] = as.data.frame(table(gbmMafs[[i]]$Variant_Classification))
  colnames(allVars_table[[i]])[2] = paste(names(gbmMafs)[i])
  names(allVars_table)[i] = names(gbmMafs)[i]
}

df_allVarsTable = allVars_table %>% reduce(full_join,by=c('Var1'))
rownames(df_allVarsTable) = df_allVarsTable$Var1
df_allVarsTable = df_allVarsTable[,c(4,2,5,3,8,6,9,7,12,11,13,10)]
temp = t(data.frame(colSums(df_allVarsTable,na.rm = TRUE)))
rownames(temp) = 'Total'
df_allVarsTable = rbind(df_allVarsTable,temp)
df_allVarsTable[is.na(df_allVarsTable)] = ''



# tables of variants per sample (FILTERED)
numRows = lapply(filtMafs,function(x) {unlist(table(x$Variant_Classification))})
numRows = max(unlist(lapply(numRows, nrow)))
filtTable = data.frame(matrix(ncol=12,nrow=numRows))

# make list of table(variant_classification) for each sample
filtTable = list()
for (i in 1:length(filtMafs)) {
  if (nrow(filtMafs[[i]]) > 0) {
    filtTable[[i]] = as.data.frame(table(filtMafs[[i]]$Variant_Classification))
    colnames(filtTable[[i]])[2] = paste(names(filtMafs)[i])
    names(filtTable)[i] = names(filtMafs)[i]
  } else {#                                                   needed because some have zero variants
    filtTable[[i]] = data.frame(matrix(nrow=1,ncol=2))
    colnames(filtTable[[i]]) = c('Var1',paste(names(filtMafs)[i]))
    filtTable[[i]]$Var1 = as.factor(filtTable[[i]]$Var1)
    filtTable[[i]][,2] = as.integer(filtTable[[i]][,2])
    names(filtTable)[i] = names(filtMafs)[i]
  }
  
}

# turn list into data frame
df_filtTable = filtTable %>% reduce(full_join,by=c('Var1'))
df_filtTable = df_filtTable[!is.na(df_filtTable$Var1),]
rownames(df_filtTable) = df_filtTable$Var1
df_filtTable = df_filtTable[,c(4,2,5,3,8,6,9,7,12,11,13,10)]
temp = t(data.frame(colSums(df_filtTable,na.rm = TRUE)))
rownames(temp) = 'Total'
df_filtTable = rbind(df_filtTable,temp)
df_filtTable[is.na(df_filtTable)] = ''



#### Set up copy number data for later ####
# read in CNVs
temp = list.files(path = '/Users/valdezkm/Documents/Cambridge/sequenza',pattern = '*segments.txt',full.names = TRUE)
allCNVs = lapply(temp, read.delim, stringsAsFactors = F)
names(allCNVs) = basename(temp)
pheno = data.frame(fileNames = basename(temp))
pheno$fileNames = gsub('.*\\+|_.*','',pheno$fileNames)
seqSampleKey = read.delim('/Users/valdezkm/Documents/Cambridge/SequencingSampleKey_McAbeeJH.txt')
pheno = merge(pheno,seqSampleKey,by.x='fileNames',by.y='File.Letter.Name',all.x=TRUE, sort=FALSE)
pheno$Label.Name = gsub('DNA|\\*|12/16|12/2016','',pheno$Label.Name)
pheno$Label.Name = gsub(' ','',pheno$Label.Name)
pheno$Label.Name = gsub('_\\(2\\)','',pheno$Label.Name)

# make sure all above worked right:
pheno$fileNames == gsub('.*\\+|_.*','',names(allCNVs))

# rename maf names
names(allCNVs) = pheno$Label.Name

# remove any NA rows
allCNVs = lapply(allCNVs,na.omit)

# assign copy numbers to mutations
assignCNVtoSNV = function(snv,cnv) {
  snv = snv[[1]]
  cnv = cnv[[1]]
  snv$CNt = NA
  snv$A = NA
  snv$B = NA
  for (i in 1:nrow(snv)) {
    for (j in 1:nrow(cnv)) {
      if ((snv$Chromosome[i] == cnv$chromosome[j]) & (snv$Start_Position[i] >= cnv$start.pos[j]) & (snv$End_Position[i] <= cnv$end.pos[j])) {
        snv$CNt[i] = cnv$CNt[j]
        snv$A[i] = cnv$A[j]
        snv$B[i] = cnv$B[j]
      }
    }
  }
  return(snv)
}

SNV_CNV = c()
for (k in 1:length(allMafs)) {
  SNV_CNV[[k]] = assignCNVtoSNV(allMafs[k],allCNVs[k])
}
names(SNV_CNV) = pheno$Label.Name


justCNs = c()
for (i in 1:length(SNV_CNV)) {
  justCNs[[i]] = subset(SNV_CNV[[i]],select=c(Hugo_Symbol,Chromosome,Start_Position,CNt))
  colnames(justCNs[[i]])[4] = names(SNV_CNV[i])
}

df_justCNs = justCNs %>% reduce(full_join, by = c('Hugo_Symbol','Chromosome','Start_Position'))
mat_justCNs = as.matrix(df_justCNs[,4:ncol(df_justCNs)])






#### Read in river plot info - use clonality for nei distance ####
J3_river = read.delim('/Users/valdezkm/Documents/Cambridge/superFreq/myAnalysis/plots/J3/rivers/J3-river.csv',sep=',')
J3_river = J3_river[J3_river$clone!='germline',]
#
J7_river = read.delim('/Users/valdezkm/Documents/Cambridge/superFreq/myAnalysis/plots/J7/rivers/J7-river.csv',sep=',')
J7_river = J7_river[J7_river$clone!='germline',]
#
J14_river = read.delim('/Users/valdezkm/Documents/Cambridge/superFreq/myAnalysis/plots/J14/rivers/J14-river.csv',sep=',')
J14_river = J14_river[J14_river$clone!='germline',]

riverList = list(J3_river=J3_river,J7_river=J7_river,J14_river=J14_river)

for (i in 1:length(riverList)) {
  riverList[[i]]$name = gsub(' \\(.*','',riverList[[i]]$name)
  riverList[[i]]$genePos = paste(riverList[[i]]$name,riverList[[i]]$start,riverList[[i]]$end)
  colnames(riverList[[i]]) = gsub('clonality\\.','',colnames(riverList[[i]]))
  riverList[[i]] = riverList[[i]][,c(22,6:9)]
  riverList[[i]] = melt(riverList[[i]])
  colnames(riverList[[i]])[2] = 'Sample'
  riverList[[i]]$Pop = 'none'
  riverList[[i]]$Ploidy = 2
  riverList[[i]]$Format = 'freq'
  riverList[[i]]$mash = paste(riverList[[i]]$genePos,riverList[[i]]$Sample)
  dupl = duplicated(riverList[[i]]$mash) | duplicated(riverList[[i]]$mash,fromLast = TRUE)   # remove duplicate variants
  riverList[[i]] = riverList[[i]][!dupl,]
  riverList[[i]] = dcast(riverList[[i]],Sample+Pop+Ploidy+Format~genePos,value.var = 'value')
}







#### Calculate nei genetic distance with first level filtered data (germline variants removed) ####
library(StAMPP)
library(reshape2)


stamps = list()
for (i in 1:length(allMafs)) {
  stamps[[i]] = allMafs[[i]]
  stamps[[i]]$genePos = paste(stamps[[i]]$Hugo_Symbol,stamps[[i]]$Start_Position,stamps[[i]]$End_Position)
  stamps[[i]] = subset(stamps[[i]],select=c(genePos,vaf))
  stamps[[i]]$Sample = names(allMafs)[i]
  stamps[[i]]$Pop = 'none'
  stamps[[i]]$Format = 'freq'
  stamps[[i]]$Ploidy = 2
  names(stamps)[i] = names(allMafs)[i]
}

# J3 Nei 
stamps_J3 = stamps[grep('J3',names(stamps))]
df_stamps_J3 = stamps_J3 %>% reduce(full_join, by = c('genePos','Sample','Pop','Format','Ploidy','vaf'))
df_stamps_J3 = dcast(df_stamps_J3,Sample+Pop+Ploidy+Format~genePos,value.var = 'vaf')
df_stamps_J3 = t(na.omit(t(df_stamps_J3)))
df_stamps_J3 = as.data.frame(df_stamps_J3)

df_stamps_J3 = stamppConvert(df_stamps_J3,type='r')

nei_J3 = data.frame(stamppNeisD(df_stamps_J3,FALSE))
colnames(nei_J3) = rownames(nei_J3)
pheatmap(nei_J3,display_numbers = T,fontsize = 14)


# J7 Nei
stamps_J7 = stamps[grep('J7',names(stamps))]
df_stamps_J7 = stamps_J7 %>% reduce(full_join, by = c('genePos','Sample','Pop','Format','Ploidy','vaf'))
df_stamps_J7 = dcast(df_stamps_J7,Sample+Pop+Ploidy+Format~genePos,value.var = 'vaf')
df_stamps_J7 = t(na.omit(t(df_stamps_J7)))
df_stamps_J7 = as.data.frame(df_stamps_J7)

df_stamps_J7 = stamppConvert(df_stamps_J7,type='r')

nei_J7 = data.frame(stamppNeisD(df_stamps_J7,FALSE))
colnames(nei_J7) = rownames(nei_J7)
pheatmap(nei_J7,display_numbers = T,fontsize = 14)

# J14 Nei
stamps_J14 = stamps[grep('J14',names(stamps))]
df_stamps_J14 = stamps_J14 %>% reduce(full_join, by = c('genePos','Sample','Pop','Format','Ploidy','vaf'))
df_stamps_J14 = dcast(df_stamps_J14,Sample+Pop+Ploidy+Format~genePos,value.var = 'vaf')
df_stamps_J14 = t(na.omit(t(df_stamps_J14)))
df_stamps_J14 = as.data.frame(df_stamps_J14)

df_stamps_J14 = stamppConvert(df_stamps_J14,type='r')

nei_J14 = data.frame(stamppNeisD(df_stamps_J14,FALSE))
colnames(nei_J14) = rownames(nei_J14)
pheatmap(nei_J14,display_numbers = T,fontsize = 14)







#### Calculate nei genetic distance with second level filtered data (germline and non-protein-altering variants removed) ####
library(StAMPP)
library(reshape2)


stamps = list()
for (i in 1:length(filtMafs)) {
  stamps[[i]] = filtMafs[[i]]
  stamps[[i]]$genePos = paste(stamps[[i]]$Hugo_Symbol,stamps[[i]]$Start_Position,stamps[[i]]$End_Position)
  stamps[[i]] = subset(stamps[[i]],select=c(genePos,vaf))
  stamps[[i]]$Sample = names(filtMafs)[i]
  stamps[[i]]$Pop = 'none'
  stamps[[i]]$Format = 'freq'
  stamps[[i]]$Ploidy = 2
  names(stamps)[i] = names(filtMafs)[i]
}

# J3 Nei 
stamps_J3 = stamps[grep('J3',names(stamps))]
df_stamps_J3 = stamps_J3 %>% reduce(full_join, by = c('genePos','Sample','Pop','Format','Ploidy','vaf'))
df_stamps_J3 = dcast(df_stamps_J3,Sample+Pop+Ploidy+Format~genePos,value.var = 'vaf')
df_stamps_J3 = t(na.omit(t(df_stamps_J3)))
df_stamps_J3 = as.data.frame(df_stamps_J3)

df_stamps_J3 = stamppConvert(df_stamps_J3,type='r')

nei_J3 = data.frame(stamppNeisD(df_stamps_J3,FALSE))
colnames(nei_J3) = rownames(nei_J3)
pheatmap(nei_J3,display_numbers = T,fontsize = 14)


# J7 Nei
stamps_J7 = stamps[grep('J7',names(stamps))]
df_stamps_J7 = stamps_J7 %>% reduce(full_join, by = c('genePos','Sample','Pop','Format','Ploidy','vaf'))
df_stamps_J7 = dcast(df_stamps_J7,Sample+Pop+Ploidy+Format~genePos,value.var = 'vaf')
df_stamps_J7 = t(na.omit(t(df_stamps_J7)))
df_stamps_J7 = as.data.frame(df_stamps_J7)

df_stamps_J7 = stamppConvert(df_stamps_J7,type='r')

nei_J7 = data.frame(stamppNeisD(df_stamps_J7,FALSE))
colnames(nei_J7) = rownames(nei_J7)
pheatmap(nei_J7,display_numbers = T,fontsize = 14)

# J14 Nei
stamps_J14 = stamps[grep('J14',names(stamps))]
df_stamps_J14 = stamps_J14 %>% reduce(full_join, by = c('genePos','Sample','Pop','Format','Ploidy','vaf'))
df_stamps_J14 = dcast(df_stamps_J14,Sample+Pop+Ploidy+Format~genePos,value.var = 'vaf')
df_stamps_J14 = t(na.omit(t(df_stamps_J14)))
df_stamps_J14 = as.data.frame(df_stamps_J14)

df_stamps_J14 = stamppConvert(df_stamps_J14,type='r')

nei_J14 = data.frame(stamppNeisD(df_stamps_J14,FALSE))
colnames(nei_J14) = rownames(nei_J14)
pheatmap(nei_J14,display_numbers = T,fontsize = 14)




























# take a look at superFreq CNA calls, see if they can be used for expands

library(superFreq)
library(expands)

J3_data = loadData('/Users/valdezkm/Documents/Cambridge/superFreq/myAnalysis/R/J3')
J7_data = loadData('/Users/valdezkm/Documents/Cambridge/superFreq/myAnalysis/R/J7')
J14_data = loadData('/Users/valdezkm/Documents/Cambridge/superFreq/myAnalysis/R/J14')

# # J3
# J3_T1_cellLine_CNA = J3_data$clusters$J3_T1_cells$clusters
# J3_T4_cellLine_CNA = J3_data[["clusters"]][["J3_T4_cells"]][["clusters"]]
# J3_T1_patient_CNA = J3_data$clusters$J3_T1$clusters
# J3_T4_patient_CNA = J3_data$clusters$J3_T4$clusters
# #
# # J7
# J7_T1_cellLine_CNA = J7_data$clusters$J7_T1_cells$clusters
# J7_T3_cellLine_CNA = J7_data$clusters$J7_T3_cells$clusters
# J7_T1_patient_CNA = J7_data$clusters$J7_T1$clusters
# J7_T3_patient_CNA = J7_data$clusters$J7_T3$clusters
# #
# # J14
# J14_T3_cellLine_CNA = J14_data$clusters$J14_T3_cells$clusters
# J14_T6_cellLine_CNA = J14_data$clusters$J14_T6_cells$clusters
# J14_T3_patient_CNA = J14_data$clusters$J14_T3$clusters
# J14_T6_patient_CNA = J14_data$clusters$J14_T6$clusters

# parse out CNs to use with expands
allData = list(D_S1 = J3_data$clusters$J3_T1_cells$clusters,
               G_S6 = J3_data[["clusters"]][["J3_T4_cells"]][["clusters"]],
               J_S9 = J3_data$clusters$J3_T1$clusters,
               K_S8 = J3_data$clusters$J3_T4$clusters,
               #
               # J7
               A_S3 = J7_data$clusters$J7_T1_cells$clusters,
               H_S5 = J7_data$clusters$J7_T3_cells$clusters,
               L_S7 = J7_data$clusters$J7_T1$clusters,
               M_S12 = J7_data$clusters$J7_T3$clusters,
               #
               # J14
               I_S4 = J14_data$clusters$J14_T3_cells$clusters,
               C_S2 = J14_data$clusters$J14_T6_cells$clusters,
               N_S11 = J14_data$clusters$J14_T3$clusters,
               O_S10 = J14_data$clusters$J14_T6$clusters)

for (i in 1:length(allData)) {
  data(cbs)
  allData[[i]][,20][allData[[i]][,20] == 'CL'] = ''
  allData[[i]]$CN = nchar(allData[[i]][,20])
  allData[[i]]$chr = gsub('\\..*','',row.names(allData[[i]]))
  allData[[i]]$chr = gsub('chr','',allData[[i]]$chr)
  allData[[i]] = subset(allData[[i]],select=c(chr,x1,x2,CN))
  colnames(allData[[i]]) = colnames(cbs)
  allData[[i]]$chr = gsub('chr','',allData[[i]]$chr)
  allData[[i]] = allData[[i]][allData[[i]]$chr!='X',]
  allData[[i]] = allData[[i]][allData[[i]]$chr!='Y',]
  allData[[i]] = sapply(allData[[i]],as.numeric)
  write.csv(allData[[i]],paste0('/Users/valdezkm/Documents/Cambridge/superFreq/CNs/',names(allData[i]),'.csv'),row.names = FALSE)
}


stories = J3_data$stories$stories$J3$clusters$cloneStories
all = J3_data$stories$stories$J3$all
snps = J3_data$allVariants$variants$SNPs
plotStories(stories,snps)





# Read in mafs to determine avg vaf in response to Andor github
temp = list.files(path = '/Users/valdezkm/Documents/Cambridge/mutect2/',pattern="*.maf",full.names = TRUE)
mafs = lapply(temp, read.delim)
names(mafs) = basename(temp)

mafAvgs = c()
for (i in 1:length(mafs)) {
  mafs[[i]]$vaf = mafs[[i]]$t_alt_count / mafs[[i]]$t_depth
  tmp = data.frame(mean(mafs[[i]]$vaf))
  rownames(tmp) = names(mafs[i])
  mafAvgs = rbind(mafAvgs,tmp)
}

# Read in sequenza CN to diagnose expands problem
temp = list.files(path = '/Users/valdezkm/Documents/Cambridge/sequenza',full.names = TRUE, pattern = '*segments.txt')
cns = lapply(temp, read.delim)
names(cns) = basename(temp)

medianCN = c()
for (i in 1:length(cns)) {
  tmp = data.frame(median(cns[[i]]$CNt,na.rm = TRUE))
  rownames(tmp) = names(cns[i])
  medianCN = rbind(medianCN,tmp)
}

meanCN = c()
for (i in 1:length(cns)) {
  tmp = data.frame(mean(cns[[i]]$CNt,na.rm = TRUE))
  rownames(tmp) = names(cns[i])
  meanCN = rbind(meanCN,tmp)
}





#### YAPSA for cosmic mutation analysis ####
library(YAPSA)

# load old set of signatures so can use this to reorder new set
Alex_signatures_path <- paste0("ftp://ftp.sanger.ac.uk/pub/cancer/",
                               "AlexandrovEtAl/signatures.txt")
AlexInitialArtif_sig_df <- read.csv(Alex_signatures_path,header=TRUE,sep="\t")
# reformat
Alex_rownames <- paste(AlexInitialArtif_sig_df[,1],
                       AlexInitialArtif_sig_df[,2],sep=" ")
select_ind <- grep("Signature",names(AlexInitialArtif_sig_df))
AlexInitialArtif_sig_df <- AlexInitialArtif_sig_df[,select_ind]
number_of_Alex_sigs <- dim(AlexInitialArtif_sig_df)[2]
names(AlexInitialArtif_sig_df) <- gsub("Signature\\.","A",
                                       names(AlexInitialArtif_sig_df))
rownames(AlexInitialArtif_sig_df) <- Alex_rownames
# narrow down old list from 27 to 22
AlexInitialValid_sig_df <- AlexInitialArtif_sig_df[,grep("^A[0-9]+",
                                                         names(AlexInitialArtif_sig_df))]
number_of_Alex_validated_sigs <- dim(AlexInitialValid_sig_df)[2]


# Load UPDATED set of mutational signatures
Alex_COSMIC_signatures_path <- 
  paste0("http://cancer.sanger.ac.uk/cancergenome/",
         "assets/signatures_probabilities.txt")
AlexCosmicValid_sig_df <- read.csv(Alex_COSMIC_signatures_path,
                                   header=TRUE,sep="\t")
Alex_COSMIC_rownames <- paste(AlexCosmicValid_sig_df[,1],
                              AlexCosmicValid_sig_df[,2],sep=" ")
COSMIC_select_ind <- grep("Signature",names(AlexCosmicValid_sig_df))
AlexCosmicValid_sig_df <- AlexCosmicValid_sig_df[,COSMIC_select_ind]
number_of_Alex_COSMIC_sigs <- dim(AlexCosmicValid_sig_df)[2]
names(AlexCosmicValid_sig_df) <- gsub("Signature\\.","AC",
                                      names(AlexCosmicValid_sig_df))
rownames(AlexCosmicValid_sig_df) <- Alex_COSMIC_rownames
# reorder features
COSMIC_order_ind <- match(Alex_rownames,Alex_COSMIC_rownames)
AlexCosmicValid_sig_df <- AlexCosmicValid_sig_df[COSMIC_order_ind,]


# preparation for later analysis
signature_colour_vector <- c("darkgreen","green","pink","goldenrod",
                             "lightblue","blue","orangered","yellow",
                             "orange","brown","purple","red",
                             "darkblue","magenta","maroon","yellowgreen",
                             "violet","lightgreen","sienna4","deeppink",
                             "darkorchid","seagreen","grey10","grey30",
                             "grey50","grey70","grey90")
bio_process_vector <- c("spontaneous deamination","spontaneous deamination",
                        "APOBEC","BRCA1_2","Smoking","unknown",
                        "defect DNA MMR","UV light exposure","unknown",
                        "IG hypermutation","POL E mutations","temozolomide",
                        "unknown","APOBEC","unknown","unknown","unknown",
                        "unknown","unknown","unknown","unknown","unknown",
                        "nonvalidated","nonvalidated","nonvalidated",
                        "nonvalidated","nonvalidated")
AlexInitialArtif_sigInd_df <- data.frame(sig=colnames(AlexInitialArtif_sig_df))
AlexInitialArtif_sigInd_df$index <- seq_len(dim(AlexInitialArtif_sigInd_df)[1])
AlexInitialArtif_sigInd_df$colour <- signature_colour_vector
AlexInitialArtif_sigInd_df$process <- bio_process_vector

COSMIC_signature_colour_vector <- c("green","pink","goldenrod",
                                    "lightblue","blue","orangered","yellow",
                                    "orange","brown","purple","red",
                                    "darkblue","magenta","maroon",
                                    "yellowgreen","violet","lightgreen",
                                    "sienna4","deeppink","darkorchid",
                                    "seagreen","grey","darkgrey",
                                    "black","yellow4","coral2","chocolate2",
                                    "navyblue","plum","springgreen")
COSMIC_bio_process_vector <- c("spontaneous deamination","APOBEC",
                               "defect DNA DSB repair hom. recomb.",
                               "tobacco mutatgens, benzo(a)pyrene",
                               "unknown",
                               "defect DNA MMR, found in MSI tumors",
                               "UV light exposure","unknown","POL eta and SHM",
                               "altered POL E",
                               "alkylating agents, temozolomide",
                               "unknown","APOBEC","unknown",
                               "defect DNA MMR","unknown","unknown",
                               "unknown","unknown",
                               "associated w. small indels at repeats",
                               "unknown","aristocholic acid","unknown",
                               "aflatoxin","unknown","defect DNA MMR",
                               "unknown","unknown","tobacco chewing","unknown")
AlexCosmicValid_sigInd_df <- data.frame(sig=colnames(AlexCosmicValid_sig_df))
AlexCosmicValid_sigInd_df$index <- seq_len(dim(AlexCosmicValid_sigInd_df)[1])
AlexCosmicValid_sigInd_df$colour <- COSMIC_signature_colour_vector
AlexCosmicValid_sigInd_df$process <- COSMIC_bio_process_vector


# Before running this had to edit the chr out of the bed and fasta files!
# Edit chr out of bed file
# $ sed -i -e 's/chr//g' v5.bed 
# 
# Edit >chr out of fasta file
# $ sed -i -e 's/>chr/>/g' genome.fa 

# This didn't work, wrong k-mer frequencies in ref likely caused by extra characters on each header line (>)
# $ sed -i -e '/^>/ s/ .*//' file.fasta 

# Still no, the actual issue was random letters in the genome... K,Y etc
# Manually take these out of kmer_frequencies_in_ref.csv then run again:
# library(stringr)
# kmer = read.delim('/Users/valdezkm/Documents/ccbrpipeliner_WES_genome_YAPSA_hg38/projectFolder/kmer_frequencies_in_ref_toBeEdited.csv',header = F)
# kmer = kmer[str_count(kmer$V1,'A|C|T|G')>2,]
# kmer = paste0(kmer$V1,'\t',kmer$V2)
# write.table(kmer,'/Users/valdezkm/Documents/ccbrpipeliner_WES_genome_YAPSA_hg38/projectFolder/kmer_frequencies_in_ref.csv',row.names = F,col.names = F,quote = F)
# 
# # Adjust genome to account for using WES instead of WGS, produce normalized correlation factors called correls
correls_all = run_kmer_frequency_correction(in_ref_genome_fasta = '/Users/valdezkm/Documents/ccbrpipeliner_WES_genome_YAPSA_hg38/GRCh38.d1.vd1.fa',
                              in_target_capture_bed = '/Users/valdezkm/Documents/ccbrpipeliner_WES_genome_YAPSA_hg38/V5_noUTR_hg38.bed',
                              in_word_length = 3,
                              project_folder = '/Users/valdezkm/Documents/ccbrpipeliner_WES_genome_YAPSA_hg38/projectFolder',
                              target_capture_fasta = 'WES_targetCapture.fa')
# correls = correls_all$rel_cor
# save(correls,file='/Volumes/BioDiscovery/Lab/Tofilon/YAPSA/correls.RData')
load('/Volumes/BioDiscovery/Lab/Tofilon/YAPSA/correls.RData')



####-------!!!!!~~~~~~ Begin enter data for YAPSA -------!!!!!~~~~~~####
# FILTERED DATA --- prepare data - change name here to change patient
yap_name = 'J14'
yap = filtMafs[grepl(yap_name,names(filtMafs))]
for (i in 1:length(yap)) {
  yap[[i]]$Sample = names(yap[i])
  yap[[i]] = yap[[i]][yap[[i]]$Variant_Type=='SNP',]
  yap[[i]] = yap[[i]][which(colnames(yap[[i]])%in%c('Chromosome','Start_Position','Reference_Allele','Tumor_Seq_Allele2','Sample','chrStart'))]
}

# UNFILTERED DATA
# yap_name = 'J14'
# yap = allMafs[grepl(yap_name,names(allMafs))]
# for (i in 1:length(yap)) {
#   yap[[i]]$Sample = names(yap[i])
#   yap[[i]] = yap[[i]][yap[[i]]$Variant_Type=='SNP',]
#   yap[[i]] = yap[[i]][which(colnames(yap[[i]])%in%c('Chromosome','Start_Position','Reference_Allele','Tumor_Seq_Allele2','Sample','chrStart'))]
# }
####-------!!!!!~~~~~~ End enter data for YAPSA -------!!!!!~~~~~~####


yap = yap %>% purrr::reduce(full_join,by=c('Chromosome','Start_Position','Reference_Allele','Tumor_Seq_Allele2','Sample'))
colnames(yap) = c("CHROM","POS","REF","ALT","PID")


## reformatting per tutorial
yap$change <- 
  attribute_nucleotide_exchanges(yap)
yap <- yap[order(yap$PID,yap$CHROM,yap$POS),]
yap <- annotate_intermut_dist_cohort(yap,in_PID.field="PID")
data("exchange_colour_vector")
yap$col <- exchange_colour_vector[yap$change]
yap$SUBGROUP = yap_name



# Building a mutational catalogue
library(BSgenome.Hsapiens.UCSC.hg38)   
word_length <- 3
yap_mutCat_list <- 
  create_mutation_catalogue_from_df(
    yap,
    this_seqnames.field = "CHROM", this_start.field = "POS",
    this_end.field = "POS", this_PID.field = "PID",
    this_subgroup.field = "SUBGROUP",
    this_refGenome = BSgenome.Hsapiens.UCSC.hg38,
    this_wordLength = word_length)
yap_mutCat_df = as.data.frame(yap_mutCat_list$matrix)

## correct counts to account for WES target capture
yap_mutCat_df = normalizeMotifs_otherRownames(in_matrix = yap_mutCat_df,
                                              in_norms = correls,
                                              adjust_counts = TRUE)

############### WITH SIGNATURE-SPECIFIC CUTOFF ############### 
# load signature-specific cut-offs
data(cutoffs)
specific_cutoff_vector = cutoffCosmicValid_rel_df[6,]

# LCD_complex_cutoff --- WITH signature-specific cutoffs
current_sig_df <- AlexCosmicValid_sig_df
current_sigInd_df <- AlexCosmicValid_sigInd_df
CosmicValid_cutoffSpec_LCDlist <- LCD_complex_cutoff(
  in_mutation_catalogue_df = yap_mutCat_df,
  in_signatures_df = current_sig_df,
  in_cutoff_vector = specific_cutoff_vector,
  in_sig_ind_df = current_sigInd_df)

# # Some adaptation (extracting and reformatting the information which sample belongs to which subgroup):
yap_COSMICExposures_df <-
  LCD(yap_mutCat_df,current_sig_df)
#
COSMIC_subgroups_df <-
  make_subgroups_df(yap,
                    yap_COSMICExposures_df)

# Plotting absolute exposures for visualization:
exposures_barplot(
  in_exposures_df = CosmicValid_cutoffSpec_LCDlist$exposures,
  in_signatures_ind_df = CosmicValid_cutoffSpec_LCDlist$out_sig_ind_df,
  in_subgroups_df = COSMIC_subgroups_df,
  in_sum_ind = c(1,2,3,4))    #     determine order of samples to be plotted

# And relative exposures:
exposures_barplot(
  in_exposures_df = CosmicValid_cutoffSpec_LCDlist$norm_exposures,
  in_signatures_ind_df = CosmicValid_cutoffSpec_LCDlist$out_sig_ind_df,
  in_subgroups_df = COSMIC_subgroups_df,
  in_sum_ind = c(1,2,3,4))

# cluster samples based on their signature exposures
complex_heatmap_exposures(CosmicValid_cutoffSpec_LCDlist$norm_exposures,
                          COSMIC_subgroups_df,
                          CosmicValid_cutoffSpec_LCDlist$out_sig_ind_df,
                          in_data_type="norm exposures",
                          in_subgroup_colour_column="col",
                          in_method="manhattan",
                          in_subgroup_column="subgroup")













##### Read in cnvkit output - convert log2 ratios to absolute scale #####
temp = list.files(path = '/Users/valdezkm/Documents/Cambridge/cnvkit',pattern = '*.cns',full.names = TRUE)
cns = lapply(temp, read.delim, stringsAsFactors = F)
names(cns) = gsub('_calls.cns','',basename(temp))
# add estimated cn without rounding 
for (i in 1:length(cns)) {
  cns[[i]]$cn = (2^cns[[i]]$log2)*2
  write.table(cns[[i]],paste0('/Users/valdezkm/Documents/Cambridge/cnvkit/',names(cns[i]),'_callsExpands.cns'),sep='\t',quote = F,row.names = F)
}

##### !!!!!!!!!!!!! #####
#### Expands - read in number of SPs per sample ####
temp = list.files(path = '/Users/valdezkm/Documents/Cambridge/expands/orig_output',pattern = '*annotatedSPs.csv',full.names = TRUE)
SPs = lapply(temp, read.delim, stringsAsFactors = F, sep=',')
names(SPs) = gsub('_.*','', basename(temp)) 
#
pheno = data.frame(fileNames = basename(temp))
pheno$fileNames = gsub('.*\\+|_.*','',pheno$fileNames)
seqSampleKey = read.delim('/Users/valdezkm/Documents/Cambridge/SequencingSampleKey_McAbeeJH.txt')
pheno = merge(pheno,seqSampleKey,by.x='fileNames',by.y='File.Letter.Name',all.x=TRUE, sort=FALSE)
pheno$Label.Name = gsub('DNA|\\*|12/16|12/2016','',pheno$Label.Name)
pheno$Label.Name = gsub(' ','',pheno$Label.Name)
pheno$Label.Name = gsub('_\\(2\\)','',pheno$Label.Name)
#
# make sure all above worked right:
pheno$fileNames == names(SPs)
#
# rename 
names(SPs) = pheno$Label.Name
#
numSPs = data.frame(Sample = c(names(SPs)))
for (i in 1:length(SPs)) {
  numSPs$numSPs[[i]] = nrow(unique(SPs[[i]][12]))
  numSPs$Purity_Expands[[i]] = round(max(SPs[[i]]$SP),digits = 2)
}
numSPs = numSPs[order(numSPs$Sample),]
numSPs$Cellularity_Sequenza = c(0.76,1,0.73,1,0.68,1,0.88,1,0.63,1,0.35,1)
write.csv(numSPs,'/Users/valdezkm/Documents/Cambridge/expands/numSPs.csv',row.names = F)



#### merge expands output with mafs ####
# GBM driver genes
# JoeGenes = read.delim('/Users/valdezkm/Documents/Driver_Genes/Joseph_GBM_Drivers.txt')
# JoeGenes = as.character(unique(JoeGenes$GBM.Drive.Genes))

# Campbell Cancer Driver Genes
campbell = read.delim('/Volumes/BioDiscovery/Lab/Tofilon/P06_Cambridge_Joseph/Campbell_Cancer_Genes.csv',sep=',',header=F)
campbell = unlist(campbell)
campbell = campbell[campbell!='']
campbell = droplevels(campbell)
#
temp = list.files(path = '/Users/valdezkm/Documents/Cambridge/expands/orig_output',pattern = '*dm.csv',full.names = TRUE)
expands = lapply(temp, read.delim, stringsAsFactors = F, sep=',')
names(expands) = gsub('_.*','', basename(temp)) 
#
pheno = data.frame(fileNames = basename(temp))
pheno$fileNames = gsub('.*\\+|_.*','',pheno$fileNames)
seqSampleKey = read.delim('/Users/valdezkm/Documents/Cambridge/SequencingSampleKey_McAbeeJH.txt')
pheno = merge(pheno,seqSampleKey,by.x='fileNames',by.y='File.Letter.Name',all.x=TRUE, sort=FALSE)
pheno$Label.Name = gsub('DNA|\\*|12/16|12/2016','',pheno$Label.Name)
pheno$Label.Name = gsub(' ','',pheno$Label.Name)
pheno$Label.Name = gsub('_\\(2\\)','',pheno$Label.Name)
#
# make sure all above worked right:
pheno$fileNames == names(expands)
#
# rename 
names(expands) = pheno$Label.Name
#
# sort mafs and expands lists so they match
filtMafs = filtMafs[order(names(filtMafs))]
expands = expands[order(names(expands))]
names(filtMafs) == names(expands)
# take chr off maf chromosome names

# annotate expands..... keep all filtered mafs, omit non-protein-altering and low vafs with all.x=TRUE
expands_annot = list()
for (i in 1:length(expands)) {
  filtMafs[[i]]$Chromosome = gsub('chr','',filtMafs[[i]]$Chromosome)
  expands_annot[[i]] = merge(filtMafs[[i]],expands[[i]],by.x=c('Chromosome','Start_Position'),by.y=c('chr','startpos'),all.x=TRUE)
  expands_annot[[i]]$Sample = names(expands)[i]
  names(expands_annot)[i] = names(expands)[i]
}

# Create column to identify clonal or subclonal variants....
addCloneStatus = function(maf) {
  orderSPs = data.frame(SP=unique(maf$SP),SP_rank=rank(unique(-maf$SP)))            #create ranking of subpopulations
  orderSPs$SP_rank[is.na(orderSPs$SP)] = NA
  orderSPs$CloneStatus[!is.na(orderSPs$SP)] = 'Subclonal'
  orderSPs$CloneStatus[orderSPs$SP_rank==1] = 'Clonal'
  #
  maf = merge(orderSPs,maf,by=c('SP'),all.y=TRUE)                                  #merge to add ranking to annotated expands output
  #maf$geneClone[!is.na(maf$SP)] = paste(maf$Hugo_Symbol[!is.na(maf$SP)],'-',maf$CloneStatus[!is.na(maf$SP)])
  maf$geneClone = paste(maf$Hugo_Symbol,'-',maf$CloneStatus)
  maf = subset(maf,select=c(SP,SP_rank,CloneStatus,geneClone,Hugo_Symbol,HGVSc,Variant_Type,Variant_Classification,Start_Position,Chromosome,vaf,Gene,dbSNP_RS,Sample))
  maf
}
# for loop to apply above function to samples
for (i in 1:length(expands_annot)){
  expands_annot[[i]] = addCloneStatus(expands_annot[[i]])
}
#
# Subset patients
J3_expands = expands_annot[grepl('J3',names(expands_annot))]
J7_expands = expands_annot[grepl('J7',names(expands_annot))]
J14_expands = expands_annot[grepl('J14',names(expands_annot))]
#
# Subset cell lines
J3_CL = J3_expands[grepl('cells',names(J3_expands))]
J7_CL = J7_expands[grepl('cells',names(J7_expands))]
J14_CL = J14_expands[grepl('cells',names(J14_expands))]

# 
library(purrr)
library(reshape2)
library(pheatmap)
# turn lists into data frames 
J3_exp_df = J3_expands %>% reduce(full_join)
J7_exp_df = J7_expands %>% reduce(full_join)
J14_exp_df = J14_expands %>% reduce(full_join)
# 
# cell lines only - turn lists into data frames
J3_CL_df = J3_CL %>% reduce(full_join)
J7_CL_df = J7_CL %>% reduce(full_join)
J14_CL_df = J14_CL %>% reduce(full_join)
#
# function to prepare data for heatmap - change out patients here
heat = function(maf) {
  maf = dcast(maf,geneClone+Hugo_Symbol+CloneStatus~Sample,length)
  maf[,4:ncol(maf)][maf[,4:ncol(maf)]>1] = 1
  for (i in 1:nrow(maf)) {
    if (grepl('Clonal',maf$CloneStatus[i])) {
      maf[i,4:ncol(maf)][which(maf[i,4:ncol(maf)]!=0)] = 1
    } else if (grepl('Subclonal',maf$CloneStatus[i])) {
      maf[i,4:ncol(maf)][which(maf[i,4:ncol(maf)]!=0)] = 3
    } else if (is.na(maf$CloneStatus)[i]){
      maf[i,4:ncol(maf)][which(maf[i,4:ncol(maf)]!=0)] = NA
    }
  }
  maf = subset(maf,select=-c(geneClone,CloneStatus))
  maf = aggregate(. ~ Hugo_Symbol,data=maf,FUN=sum,na.action = na.pass)
  rownames(maf) = maf$Hugo_Symbol
  maf = subset(maf,select = -Hugo_Symbol)
  # comment out below if need df
  maf = as.matrix(maf)
}
J3_heat = heat(J3_exp_df)
J7_heat = heat(J7_exp_df)
J14_heat = heat(J14_exp_df)
# 
# cell lines only
J3_heat = heat(J3_CL_df)
J7_heat = heat(J7_CL_df)
J14_heat = heat(J14_CL_df)

# J3 Filtered Variants
setEPS()
postscript('/Users/valdezkm/Documents/Cambridge/expands/J3_heatmap_Filtered.eps', height = 20, width = 10)
pheatmap(J3_heat,na_col = 'gray',drop_levels = T,cluster_rows = FALSE, cluster_cols = FALSE, main = 'J3 - All Variants (Filtered)',cellheight = 10,cellwidth = 10,color = c('white','red','blue','springgreen4'),breaks = c(0,0.5,1.5,3.5,5))
dev.off()
#
# J3 GBM variants
J3_heat_GBM = J3_heat[rownames(J3_heat) %in% JoeGenes,]
setEPS()
postscript('/Users/valdezkm/Documents/Cambridge/expands/J3_heatmap_Filtered_GBM.eps')
pheatmap(J3_heat_GBM,na_col = 'gray',drop_levels = T,cluster_rows = FALSE, cluster_cols = FALSE, main = 'J3 - GBM Driver Genes (Filtered)',cellheight = 10,cellwidth = 10,color = c('white','red','blue','springgreen4'),breaks = c(0,0.5,1.5,3.5,5))
dev.off()
#
# J7 Filtered Variants
setEPS()
postscript('/Users/valdezkm/Documents/Cambridge/expands/J7_heatmap_Filtered.eps', height = 30, width = 10)
pheatmap(J7_heat,na_col = 'gray',drop_levels = T,cluster_rows = FALSE, cluster_cols = FALSE, main = 'J7 - All Variants (Filtered)',cellheight = 10,cellwidth = 10,color = c('white','red','blue','springgreen4'),breaks = c(0,0.5,1.5,3.5,5))
dev.off()
#
# J7 GBM variants
J7_heat_GBM = J7_heat[rownames(J7_heat) %in% JoeGenes,]
setEPS()
postscript('/Users/valdezkm/Documents/Cambridge/expands/J7_heatmap_Filtered_GBM.eps')
pheatmap(J7_heat_GBM,na_col = 'gray',drop_levels = T,cluster_rows = FALSE, cluster_cols = FALSE, main = 'J7 - GBM Driver Genes (Filtered)',cellheight = 10,cellwidth = 10,color = c('white','red','blue','springgreen4'),breaks = c(0,0.5,1.5,3.5,5))
dev.off()
#
# J14 Filtered Variants
setEPS()
postscript('/Users/valdezkm/Documents/Cambridge/expands/J14_heatmap_Filtered.eps', height = 20, width = 10)
pheatmap(J14_heat,na_col = 'gray',drop_levels = T,cluster_rows = FALSE, cluster_cols = FALSE, main = 'J14 - All Variants (Filtered)',cellheight = 10,cellwidth = 10,color = c('white','red','blue','springgreen4'),breaks = c(0,0.5,1.5,3.5,5))
dev.off()
#
# J14 GBM variants
J14_heat_GBM = J14_heat[rownames(J14_heat) %in% JoeGenes,]
setEPS()
postscript('/Users/valdezkm/Documents/Cambridge/expands/J14_heatmap_Filtered_GBM.eps')
pheatmap(J14_heat_GBM,na_col = 'gray',drop_levels = T,cluster_rows = FALSE, cluster_cols = FALSE, main = 'J14 - GBM Driver Genes (Filtered)',cellheight = 10,cellwidth = 10,color = c('white','red','blue','springgreen4'),breaks = c(0,0.5,1.5,3.5,5))
dev.off()



# Campbell - with indiv patients

# J3 GBM variants
J3_heat = J3_heat[rownames(J3_heat) %in% campbell,]
setEPS()
postscript('/Volumes/BioDiscovery/Lab/Tofilon/P06_Cambridge_Joseph/expands/J3_heatmap_campbell.eps')
pheatmap(J3_heat,na_col = 'gray',drop_levels = T,cluster_rows = FALSE, cluster_cols = FALSE, main = 'J3 - Campbell Cancer Driver Genes',cellheight = 10,cellwidth = 10,color = c('white','red','blue','springgreen4'),breaks = c(0,0.5,1.5,3.5,5))
dev.off()

# J7 GBM variants
J7_heat = J7_heat[rownames(J7_heat) %in% campbell,]
setEPS()
postscript('/Volumes/BioDiscovery/Lab/Tofilon/P06_Cambridge_Joseph/expands/J7_heatmap_campbell.eps')
pheatmap(J7_heat,na_col = 'gray',drop_levels = T,cluster_rows = FALSE, cluster_cols = FALSE, main = 'J7 - Campbell Cancer Driver Genes',cellheight = 10,cellwidth = 10,color = c('white','red','blue','springgreen4'),breaks = c(0,0.5,1.5,3.5,5))
dev.off()

# J14 GBM variants
J14_heat = J14_heat[rownames(J14_heat) %in% campbell,]
setEPS()
postscript('/Volumes/BioDiscovery/Lab/Tofilon/P06_Cambridge_Joseph/expands/J14_heatmap_campbell.eps')
pheatmap(J14_heat,na_col = 'gray',drop_levels = T,cluster_rows = FALSE, cluster_cols = FALSE, main = 'J14 - Campbell Cancer Driver Genes',cellheight = 10,cellwidth = 10,color = c('white','red','blue','springgreen4'),breaks = c(0,0.5,1.5,3.5,5))
dev.off()







# Need to fix below!

### Campbell genes with cell lines
# combine heatmaps together then subset by cambell genes
# J3_heat$rn = rownames(J3_heat)
# J7_heat$rn = rownames(J7_heat)
# J14_heat$rn = rownames(J14_heat)
# cellLines = list(J3_heat,J7_heat,J14_heat)
# cellLines = cellLines %>% reduce(full_join,by='rn')
# cellLines = cellLines[grepl('rn|cells',colnames(cellLines))]
# rownames(cellLines) = cellLines$rn
# cellLines = subset(cellLines,select=-rn)
# 
# 
# campbell_cells = cellLines[rownames(cellLines) %in% campbell,]
# setEPS()
# postscript('/Volumes/BioDiscovery/Lab/Tofilon/P06_Cambridge_Joseph/expands/heatmap_Campbell.eps')
# pheatmap(campbell_cells,na_col = 'gray',drop_levels = T,cluster_rows = FALSE, cluster_cols = FALSE, main = 'Campbell Cancer Driver Genes (Filtered)',cellheight = 10,cellwidth = 10,color = c('white','red','blue','springgreen4'),breaks = c(0,0.5,1.5,3.5,5))
# dev.off()








#### ABSOLUTE with filtered data (only protein-coding) ####
# upload segmetrics files
temp = list.files(path = '/Users/valdezkm/Documents/Cambridge/cnvkit',pattern = '*.segmetrics',full.names = TRUE)
segs = lapply(temp, read.delim, stringsAsFactors = F)
names(segs) = basename(temp)
pheno = data.frame(fileNames = basename(temp))
pheno$fileNames = gsub('.*\\+|_.*','',pheno$fileNames)
seqSampleKey = read.delim('/Users/valdezkm/Documents/Cambridge/SequencingSampleKey_McAbeeJH.txt')
pheno = merge(pheno,seqSampleKey,by.x='fileNames',by.y='File.Letter.Name',all.x=TRUE, sort=FALSE)
pheno$Label.Name = gsub('DNA|\\*|12/16|12/2016','',pheno$Label.Name)
pheno$Label.Name = gsub(' ','',pheno$Label.Name)
pheno$Label.Name = gsub('_\\(2\\)','',pheno$Label.Name)
# make sure all above worked right:
pheno$fileNames == gsub('_.*','',names(segs))
# rename maf names
names(segs) = pheno$Label.Name

# FUNCTION to format segmented copy ratios data file
formatSeg = function(input) {
  CN = subset(input, select=(c(chromosome,start,end,probes,mean)))
  colnames(CN) = c('Chromosome','Start','End','Num_Probes','Segment_Mean')
  CN = CN[!CN$Chromosome=='MT',]
  CN$Chromosome = gsub('chr','',CN$Chromosome)
  CN
}
for (i in 1:length(segs)) {
  segs[[i]] = formatSeg(segs[[i]])
  write.table(segs[[i]],paste0('/Users/valdezkm/Documents/Cambridge/ABSOLUTE/CN/',names(segs)[i],'.abs.cn.txt'),sep='\t',quote=F,row.names=F)
}

#FUNCTION to format maf file
formatMAF = function(input) {
  subMAF = subset(input, select=c(t_ref_count,t_alt_count,dbSNP_Val_Status,Start_Position,End_Position,Tumor_Sample_Barcode,Hugo_Symbol,Chromosome))
  colnames(subMAF)[4] = 'Start_position'
  colnames(subMAF)[5] = 'End_position'
  subMAF = subMAF[!(subMAF$t_alt_count==0 & subMAF$t_ref_count==0),]      #if both are zero it breaks ABSOLUTE
  subMAF$Chromosome = gsub('chr','',subMAF$Chromosome)
  subMAF
}
abs_mafs = list()
for (i in 1:length(filtMafs)) {
  abs_mafs[[i]] = formatMAF(filtMafs[[i]])
  names(abs_mafs)[i] = names(filtMafs)[i]
  write.table(abs_mafs[[i]],paste0('/Users/valdezkm/Documents/Cambridge/ABSOLUTE/MAF/',names(abs_mafs)[i],'.abs.maf.txt'),sep='\t',quote=F,row.names=F,na='')
}
##### Read in ABSOLUTE output - Filtered (protein-coding only) ######
temp = list.files(path = '/Users/valdezkm/Documents/Cambridge/ABSOLUTE/Final_Review',pattern = '_ABS_MAF.txt',full.names = TRUE)
abs = lapply(temp, read.delim, stringsAsFactors = F)
names(abs) = gsub('_ABS_MAF.txt','',basename(temp))  

#### Nei genetic distance using CCFs from ABSOLUTE ####
library(StAMPP)
library(reshape2)
library(purrr)
library(dplyr)
library(pheatmap)
for (i in 1:length(abs)) {
  abs[[i]]$genePos = paste(abs[[i]]$Hugo_Symbol,abs[[i]]$Start_position,abs[[i]]$End_position)
  abs[[i]] = subset(abs[[i]],select=c(genePos,cancer_cell_frac))
  abs[[i]]$Sample = names(abs)[i]
  abs[[i]]$Pop = 'none'
  abs[[i]]$Format = 'freq'
  abs[[i]]$Ploidy = 2
}

# J3 Nei 
abs_J3 = abs[grep('J3',names(abs))]
df_abs_J3 = abs_J3 %>% reduce(full_join, by = c('genePos','Sample','Pop','Format','Ploidy','cancer_cell_frac'))
df_abs_J3 = dcast(df_abs_J3,Sample+Pop+Ploidy+Format~genePos,value.var = 'cancer_cell_frac')
df_abs_J3 = as.data.frame(df_abs_J3)
#
df_abs_J3 = stamppConvert(df_abs_J3,type='r')

nei_J3 = data.frame(stamppNeisD(df_abs_J3,FALSE))
colnames(nei_J3) = rownames(nei_J3)
pheatmap(nei_J3,display_numbers = T,fontsize = 14)
#
# J7 Nei 
abs_J7 = abs[grep('J7',names(abs))]
df_abs_J7 = abs_J7 %>% reduce(full_join, by = c('genePos','Sample','Pop','Format','Ploidy','cancer_cell_frac'))
df_abs_J7 = dcast(df_abs_J7,Sample+Pop+Ploidy+Format~genePos,value.var = 'cancer_cell_frac')
df_abs_J7 = as.data.frame(df_abs_J7)
#
df_abs_J7 = stamppConvert(df_abs_J7,type='r')

nei_J7 = data.frame(stamppNeisD(df_abs_J7,FALSE))
colnames(nei_J7) = rownames(nei_J7)
pheatmap(nei_J7,display_numbers = T,fontsize = 14)
#
# J14 Nei 
abs_J14 = abs[grep('J14',names(abs))]
df_abs_J14 = abs_J14 %>% reduce(full_join, by = c('genePos','Sample','Pop','Format','Ploidy','cancer_cell_frac'))
df_abs_J14 = dcast(df_abs_J14,Sample+Pop+Ploidy+Format~genePos,value.var = 'cancer_cell_frac')
df_abs_J14 = as.data.frame(df_abs_J14)
#
df_abs_J14 = stamppConvert(df_abs_J14,type='r')

nei_J14 = data.frame(stamppNeisD(df_abs_J14,FALSE))
colnames(nei_J14) = rownames(nei_J14)
pheatmap(nei_J14,display_numbers = T,fontsize = 14)





#### ABSOLUTE unfiltered ####
# upload segmetrics files
temp = list.files(path = '/Users/valdezkm/Documents/Cambridge/cnvkit',pattern = '*.segmetrics',full.names = TRUE)
segs = lapply(temp, read.delim, stringsAsFactors = F)
names(segs) = basename(temp)
pheno = data.frame(fileNames = basename(temp))
pheno$fileNames = gsub('.*\\+|_.*','',pheno$fileNames)
seqSampleKey = read.delim('/Users/valdezkm/Documents/Cambridge/SequencingSampleKey_McAbeeJH.txt')
pheno = merge(pheno,seqSampleKey,by.x='fileNames',by.y='File.Letter.Name',all.x=TRUE, sort=FALSE)
pheno$Label.Name = gsub('DNA|\\*|12/16|12/2016','',pheno$Label.Name)
pheno$Label.Name = gsub(' ','',pheno$Label.Name)
pheno$Label.Name = gsub('_\\(2\\)','',pheno$Label.Name)
# make sure all above worked right:
pheno$fileNames == gsub('_.*','',names(segs))
# rename maf names
names(segs) = pheno$Label.Name

# FUNCTION to format segmented copy ratios data file
formatSeg = function(input) {
  CN = subset(input, select=(c(chromosome,start,end,probes,mean)))
  colnames(CN) = c('Chromosome','Start','End','Num_Probes','Segment_Mean')
  CN = CN[!CN$Chromosome=='MT',]
  CN$Chromosome = gsub('chr','',CN$Chromosome)
  CN
}
for (i in 1:length(segs)) {
  segs[[i]] = formatSeg(segs[[i]])
  write.table(segs[[i]],paste0('/Users/valdezkm/Documents/Cambridge/ABSOLUTE_unfiltered/CN/',names(segs)[i],'.abs.cn.txt'),sep='\t',quote=F,row.names=F)
}

#FUNCTION to format maf file
formatMAF = function(input) {
  subMAF = subset(input, select=c(t_ref_count,t_alt_count,dbSNP_Val_Status,Start_Position,End_Position,Tumor_Sample_Barcode,Hugo_Symbol,Chromosome))
  colnames(subMAF)[4] = 'Start_position'
  colnames(subMAF)[5] = 'End_position'
  subMAF = subMAF[!(subMAF$t_alt_count==0 & subMAF$t_ref_count==0),]      #if both are zero it breaks ABSOLUTE
  subMAF$Chromosome = gsub('chr','',subMAF$Chromosome)
  subMAF
}
abs_mafs = list()
for (i in 1:length(allMafs)) {
  abs_mafs[[i]] = formatMAF(allMafs[[i]])
  names(abs_mafs)[i] = names(allMafs)[i]
  write.table(abs_mafs[[i]],paste0('/Users/valdezkm/Documents/Cambridge/ABSOLUTE_unfiltered/MAF/',names(abs_mafs)[i],'.abs.maf.txt'),sep='\t',quote=F,row.names=F,na='')
}


##### Read in ABSOLUTE output - unfiltered ######
temp = list.files(path = '/Users/valdezkm/Documents/Cambridge/ABSOLUTE_unfiltered/Final_Review',pattern = '_ABS_MAF.txt',full.names = TRUE)
abs = lapply(temp, read.delim, stringsAsFactors = F)
names(abs) = gsub('_ABS_MAF.txt','',basename(temp))  

#### Nei genetic distance using CCFs from ABSOLUTE ####
library(StAMPP)
library(reshape2)
library(purrr)
library(dplyr)
library(pheatmap)
for (i in 1:length(abs)) {
  abs[[i]]$genePos = paste(abs[[i]]$Hugo_Symbol,abs[[i]]$Start_position,abs[[i]]$End_position)
  abs[[i]] = subset(abs[[i]],select=c(genePos,cancer_cell_frac))
  abs[[i]]$Sample = names(abs)[i]
  abs[[i]]$Pop = 'none'
  abs[[i]]$Format = 'freq'
  abs[[i]]$Ploidy = 2
}

# J3 Nei 
abs_J3 = abs[grep('J3',names(abs))]
df_abs_J3 = abs_J3 %>% reduce(full_join, by = c('genePos','Sample','Pop','Format','Ploidy','cancer_cell_frac'))
df_abs_J3 = dcast(df_abs_J3,Sample+Pop+Ploidy+Format~genePos,value.var = 'cancer_cell_frac')
df_abs_J3 = as.data.frame(df_abs_J3)
#
df_abs_J3 = stamppConvert(df_abs_J3,type='r')

nei_J3 = data.frame(stamppNeisD(df_abs_J3,FALSE))
colnames(nei_J3) = rownames(nei_J3)
pheatmap(nei_J3,display_numbers = T,fontsize = 14)
#
# J7 Nei 
abs_J7 = abs[grep('J7',names(abs))]
df_abs_J7 = abs_J7 %>% reduce(full_join, by = c('genePos','Sample','Pop','Format','Ploidy','cancer_cell_frac'))
df_abs_J7 = dcast(df_abs_J7,Sample+Pop+Ploidy+Format~genePos,value.var = 'cancer_cell_frac')
df_abs_J7 = as.data.frame(df_abs_J7)
#
df_abs_J7 = stamppConvert(df_abs_J7,type='r')

nei_J7 = data.frame(stamppNeisD(df_abs_J7,FALSE))
colnames(nei_J7) = rownames(nei_J7)
pheatmap(nei_J7,display_numbers = T,fontsize = 14)
#
# J14 Nei 
abs_J14 = abs[grep('J14',names(abs))]
df_abs_J14 = abs_J14 %>% reduce(full_join, by = c('genePos','Sample','Pop','Format','Ploidy','cancer_cell_frac'))
df_abs_J14 = dcast(df_abs_J14,Sample+Pop+Ploidy+Format~genePos,value.var = 'cancer_cell_frac')
df_abs_J14 = as.data.frame(df_abs_J14)
#
df_abs_J14 = stamppConvert(df_abs_J14,type='r')

nei_J14 = data.frame(stamppNeisD(df_abs_J14,FALSE))
colnames(nei_J14) = rownames(nei_J14)
pheatmap(nei_J14,display_numbers = T,fontsize = 14)



#### Copy number correlations ####
# read in CNVs
temp = list.files(path = '/Volumes/BioDiscovery/Lab/Tofilon/P06_Cambridge_Joseph/sequenza',pattern = '*segments.txt',full.names = TRUE)
allCNVs = lapply(temp, read.delim, stringsAsFactors = F)
names(allCNVs) = basename(temp)
pheno = data.frame(fileNames = basename(temp))
pheno$fileNames = gsub('.*\\+|_.*','',pheno$fileNames)
seqSampleKey = read.delim('/Users/valdezkm/Documents/Cambridge/SequencingSampleKey_McAbeeJH.txt')
pheno = merge(pheno,seqSampleKey,by.x='fileNames',by.y='File.Letter.Name',all.x=TRUE, sort=FALSE)
pheno$Label.Name = gsub('DNA|\\*|12/16|12/2016','',pheno$Label.Name)
pheno$Label.Name = gsub(' ','',pheno$Label.Name)
pheno$Label.Name = gsub('_\\(2\\)','',pheno$Label.Name)

# make sure all above worked right:
pheno$fileNames == gsub('.*\\+|_.*','',names(allCNVs))

# rename maf names
names(allCNVs) = pheno$Label.Name

# remove any NA rows
allCNVs = lapply(allCNVs,na.omit)

# assign copy numbers to mutations
assignCNVtoSNV = function(snv,cnv) {
  snv = snv[[1]]
  cnv = cnv[[1]]
  snv$CNt = NA
  snv$A = NA
  snv$B = NA
  for (i in 1:nrow(snv)) {
    for (j in 1:nrow(cnv)) {
      if ((snv$Chromosome[i] == cnv$chromosome[j]) & (snv$Start_Position[i] >= cnv$start.pos[j]) & (snv$End_Position[i] <= cnv$end.pos[j])) {
        snv$CNt[i] = cnv$CNt[j]
        snv$A[i] = cnv$A[j]
        snv$B[i] = cnv$B[j]
      }
    }
  }
  return(snv)
}

#### CN all variants ####
SNV_CNV = c()
for (k in 1:length(allCNVs)) {
  SNV_CNV[[k]] = assignCNVtoSNV(allMafs[k],allCNVs[k])
}
names(SNV_CNV) = pheno$Label.Name

justCNs = c()
for (i in 1:length(SNV_CNV)) {
  justCNs[[i]] = subset(SNV_CNV[[i]],select=c(Hugo_Symbol,Chromosome,Start_Position,CNt))
  colnames(justCNs[[i]])[4] = names(SNV_CNV[i])
}

# df_justCNs = justCNs %>% reduce(full_join, by = c('Hugo_Symbol','Chromosome','Start_Position'))
# mat_justCNs = as.matrix(df_justCNs[,4:ncol(df_justCNs)])

#J3
J3_CN = justCNs[grepl('J3',names(SNV_CNV))]
J3_CN = J3_CN %>% reduce(full_join, by = c('Hugo_Symbol','Chromosome','Start_Position'))
J3_CN = as.matrix(J3_CN[,4:ncol(J3_CN)])
cors = rcorr(J3_CN, type = 'pearson')
cors_r = cors$r
pheatmap(cors_r,display_numbers = T,fontsize = 14)

#J7
J7_CN = justCNs[grepl('J7',names(SNV_CNV))]
J7_CN = J7_CN %>% reduce(full_join, by = c('Hugo_Symbol','Chromosome','Start_Position'))
J7_CN = as.matrix(J7_CN[,4:ncol(J7_CN)])
cors = rcorr(J7_CN, type = 'pearson')
cors_r = cors$r
pheatmap(cors_r,display_numbers = T,fontsize = 14)

#J14
J14_CN = justCNs[grepl('J14',names(SNV_CNV))]
J14_CN = J14_CN %>% reduce(full_join, by = c('Hugo_Symbol','Chromosome','Start_Position'))
J14_CN = as.matrix(J14_CN[,4:ncol(J14_CN)])
cors = rcorr(J14_CN, type = 'pearson')
cors_r = cors$r
pheatmap(cors_r,display_numbers = T,fontsize = 14)







#### CN filtered variants ####
SNV_CNV = c()
for (k in 1:length(allCNVs)) {
  SNV_CNV[[k]] = assignCNVtoSNV(filtMafs[k],allCNVs[k])
}
names(SNV_CNV) = pheno$Label.Name

justCNs = c()
for (i in 1:length(SNV_CNV)) {
  justCNs[[i]] = subset(SNV_CNV[[i]],select=c(Hugo_Symbol,Chromosome,Start_Position,CNt))
  colnames(justCNs[[i]])[4] = names(SNV_CNV[i])
}

df_justCNs = justCNs %>% reduce(full_join, by = c('Hugo_Symbol','Chromosome','Start_Position'))
mat_justCNs = as.matrix(df_justCNs[,4:ncol(df_justCNs)])

#J3
J3_CN = justCNs[grepl('J3',names(SNV_CNV))]
J3_CN = J3_CN %>% reduce(full_join, by = c('Hugo_Symbol','Chromosome','Start_Position'))
J3_CN = as.matrix(J3_CN[,4:ncol(J3_CN)])
cors = rcorr(J3_CN, type = 'pearson')
cors_r = cors$r
pheatmap(cors_r,display_numbers = T,fontsize = 14)

#J7
J7_CN = justCNs[grepl('J7',names(SNV_CNV))]
J7_CN = J7_CN %>% reduce(full_join, by = c('Hugo_Symbol','Chromosome','Start_Position'))
J7_CN = as.matrix(J7_CN[,4:ncol(J7_CN)])
cors = rcorr(J7_CN, type = 'pearson')
cors_r = cors$r
pheatmap(cors_r,display_numbers = T,fontsize = 14)

#J14
J14_CN = justCNs[grepl('J14',names(SNV_CNV))]
J14_CN = J14_CN %>% reduce(full_join, by = c('Hugo_Symbol','Chromosome','Start_Position'))
J14_CN = as.matrix(J14_CN[,4:ncol(J14_CN)])
cors = rcorr(J14_CN, type = 'pearson')
cors_r = cors$r
pheatmap(cors_r,display_numbers = T,fontsize = 14)




#### Joseph thesis defense question ####
# Why is J3T4 superFreq copy number profile so different than the others? #
# check superFreq CN output #
library(superFreq)
library(reshape2)
library(dplyr)
library(purrr)
data = loadData('/Volumes/BioDiscovery/Lab/Tofilon/P06_Cambridge_Joseph/superFreq/myAnalysis/R/J3')
clusters = data$clusters
clusters = clusters[-5]   # delete norm

calls = list()
for (i in 1:length(clusters)) {
  calls[[i]] = clusters[[i]]$clusters[,'call',drop=FALSE]
  calls[[i]]$chr = row.names(calls[[i]])
  colnames(calls[[i]])[1] = names(clusters)[i]
  names(calls)[i] = names(clusters)[i]
}
calls = calls %>% reduce(full_join, by = 'chr')
calls = calls[,c(2,5,3,4,1)]

write.csv(calls, '/Volumes/BioDiscovery/Lab/Tofilon/P06_Cambridge_Joseph/superFreq_CNcalls.csv',row.names = F, na='')


#### Generate test IPA files to see if ACMG or Variant gain/loss can produce zscores in IPA #####
test = read.delim('/Users/valdezkm/Documents/GBM_cellLines_MicroArray/IPA/ipa_CTL_1321N1_Norm.txt')
test$test = runif(nrow(test),min = -2,max = 2)
write.table(test,'~/Desktop/testIPA.txt',sep='\t',row.names = F,quote = F)

test2 = read.delim('/Volumes/BioDiscovery/Lab/Tofilon/P06_Cambridge_Joseph/mutect2/P_S15+D_S1.maf')
test2 = test2[,1,drop=FALSE]
test2$test = round(runif(nrow(test2),min = -2,max = 2))
write.table(test2,'~/Desktop/testIPA.txt',sep='\t',row.names = F,quote = F)






#### Tumor only (TO) cell line (CL) analysis ####
#### Read in TO CL maf files ####
temp = list.files(path = '/Volumes/BioDiscovery/Lab/Tofilon/P06_Cambridge_Joseph/cellLines_TO/mutect2',pattern = '*.maf',full.names = TRUE)
filtTOmafs = lapply(temp, read.delim, stringsAsFactors = F)
names(filtTOmafs) = basename(temp)

for (i in 1:length(filtTOmafs)) {
  filtTOmafs[[i]]$vaf = filtTOmafs[[i]]$t_alt_count / filtTOmafs[[i]]$t_depth
  filtTOmafs[[i]] = filtTOmafs[[i]][filtTOmafs[[i]]$vaf > 0.05,]
}
#
# filter out non-protein-altering variants
variantsToKeep = c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site","Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del","In_Frame_Ins", "Missense_Mutation")
for (i in 1:length(filtTOmafs)) {
  filtTOmafs[[i]] = filtTOmafs[[i]][filtTOmafs[[i]]$Variant_Classification %in% variantsToKeep,]
}

# Write out filtered TO CL maf files for input into expands (too many variants for expands to work properly) 
# for (i in 1:length(filtTOmafs)) {
#   write.table(filtTOmafs[[i]],paste0('/Volumes/BioDiscovery/Lab/Tofilon/P06_Cambridge_Joseph/cellLines_TO/mutect2/FilteredForExpands',gsub('.maf','',names(filtTOmafs[i]),fixed=TRUE),'_filtered.maf'),quote = F,sep = '\t',row.names = F)
# }




#### CL TO: Read in expands output ####
temp = list.files(path = '/Volumes/BioDiscovery/Lab/Tofilon/P06_Cambridge_Joseph/cellLines_TO/expands_originalOutput',pattern = '*dm.csv',full.names = TRUE)
expands = lapply(temp, read.delim, stringsAsFactors = F, sep=',')
names(expands) = gsub('_.*','', basename(temp)) 
#
pheno = data.frame(fileNames = basename(temp))
pheno$fileNames = gsub('.*\\+|_.*','',pheno$fileNames)
seqSampleKey = read.delim('/Users/valdezkm/Documents/Cambridge/SequencingSampleKey_McAbeeJH.txt')
pheno = merge(pheno,seqSampleKey,by.x='fileNames',by.y='File.Letter.Name',all.x=TRUE, sort=FALSE)
pheno$Label.Name = gsub('DNA|\\*|12/16|12/2016','',pheno$Label.Name)
pheno$Label.Name = gsub(' ','',pheno$Label.Name)
pheno$Label.Name = gsub('_\\(2\\)','',pheno$Label.Name)
#
# make sure all above worked right:
pheno$fileNames == names(expands)
# rename 
names(expands) = pheno$Label.Name
#
# rename filtTOmafs
names(filtTOmafs) = gsub('_.*','', names(filtTOmafs)) 
pheno$fileNames == names(filtTOmafs)
names(filtTOmafs) = pheno$Label.Name


# Campbell Cancer Driver Genes
campbell = read.delim('/Volumes/BioDiscovery/Lab/Tofilon/P06_Cambridge_Joseph/Campbell_Cancer_Genes.csv',sep=',',header=F)
campbell = unlist(campbell)
campbell = campbell[campbell!='']
campbell = droplevels(campbell)

#
# make sure filt and expands lists match
names(filtTOmafs) == names(expands)
#
# annotate expands..... if all variants run through exands, keep all filtered mafs, omit non-protein-altering and low vafs with all.x=TRUE
expands_annot = list()
for (i in 1:length(expands)) {
  filtTOmafs[[i]]$Chromosome = gsub('chr','',filtTOmafs[[i]]$Chromosome)
  expands_annot[[i]] = merge(filtTOmafs[[i]],expands[[i]],by.x=c('Chromosome','Start_Position'),by.y=c('chr','startpos'),all.x=TRUE)
  expands_annot[[i]]$Sample = names(expands)[i]
  names(expands_annot)[i] = names(expands)[i]
}
# Function to create column to identify clonal or subclonal variants....
addCloneStatus = function(maf) {
  orderSPs = data.frame(SP=unique(maf$SP),SP_rank=rank(unique(-maf$SP)))            #create ranking of subpopulations
  orderSPs$SP_rank[is.na(orderSPs$SP)] = NA
  orderSPs$CloneStatus[!is.na(orderSPs$SP)] = 'Subclonal'
  orderSPs$CloneStatus[orderSPs$SP_rank==1] = 'Clonal'
  #
  maf = merge(orderSPs,maf,by=c('SP'),all.y=TRUE)                                  #merge to add ranking to annotated expands output
  #maf$geneClone[!is.na(maf$SP)] = paste(maf$Hugo_Symbol[!is.na(maf$SP)],'-',maf$CloneStatus[!is.na(maf$SP)])
  maf$geneClone = paste(maf$Hugo_Symbol,'-',maf$CloneStatus)
  maf = subset(maf,select=c(SP,SP_rank,CloneStatus,geneClone,Hugo_Symbol,HGVSc,Variant_Type,Variant_Classification,Start_Position,Chromosome,vaf,Gene,dbSNP_RS,Sample))
  maf
}
# for loop to apply above function to samples
for (i in 1:length(expands_annot)){
  expands_annot[[i]] = addCloneStatus(expands_annot[[i]])
}
# turn list into data frame
J3 = expands_annot[grepl('J3',names(expands_annot))] %>% reduce(full_join)
J7 = expands_annot[grepl('J7',names(expands_annot))] %>% reduce(full_join)
J14 = expands_annot[grepl('J14',names(expands_annot))] %>% reduce(full_join)
#
# function to prepare data for heatmap 
heat = function(maf) {
  maf = dcast(maf,geneClone+Hugo_Symbol+CloneStatus~Sample,length)
  maf[,4:ncol(maf)][maf[,4:ncol(maf)]>1] = 1
  for (i in 1:nrow(maf)) {
    if (grepl('Clonal',maf$CloneStatus[i])) {
      maf[i,4:ncol(maf)][which(maf[i,4:ncol(maf)]!=0)] = 1
    } else if (grepl('Subclonal',maf$CloneStatus[i])) {
      maf[i,4:ncol(maf)][which(maf[i,4:ncol(maf)]!=0)] = 3
    } else if (is.na(maf$CloneStatus)[i]){
      maf[i,4:ncol(maf)][which(maf[i,4:ncol(maf)]!=0)] = NA
    }
  }
  maf = subset(maf,select=-c(geneClone,CloneStatus))
  maf = aggregate(. ~ Hugo_Symbol,data=maf,FUN=sum,na.action = na.pass)
  rownames(maf) = maf$Hugo_Symbol
  maf = subset(maf,select = -Hugo_Symbol)
  maf = as.matrix(maf)
}
J3 = heat(J3)
J7 = heat(J7)
J14 = heat(J14)
# Filtered Variants - printed but too big for ppt
# setEPS()
# postscript('/Volumes/BioDiscovery/Lab/Tofilon/P07_A25_J14/A25/expands_output/filteredVariants_heatmap.eps', height = 110, width = 10)
# pheatmap(expands_annot,na_col = 'gray',drop_levels = T,cluster_rows = FALSE, cluster_cols = FALSE, main = 'All Variants (Filtered)',cellheight = 10,cellwidth = 10,color = c('white','red','blue','springgreen4'),breaks = c(0,0.5,1.5,3.5,5))
# dev.off()
#
# GBM variants
# heat_GBM = expands_annot[rownames(expands_annot) %in% JoeGenes,]
# setEPS()
# postscript('/Volumes/BioDiscovery/Lab/Tofilon/P07_A25_J14/A25/expands_output/GBMvariants_heatmap.eps')
# pheatmap(heat_GBM,na_col = 'gray',drop_levels = T,cluster_rows = FALSE, cluster_cols = FALSE, main = 'GBM Driver Genes (Filtered)',cellheight = 10,cellwidth = 10,color = c('white','red','blue','springgreen4'),breaks = c(0,0.5,1.5,3.5,5))
# dev.off()
#
# Campbell variants
# J3
heat_GBM = J3[rownames(J3) %in% campbell,]
setEPS()
postscript('/Volumes/BioDiscovery/Lab/Tofilon/P06_Cambridge_Joseph/cellLines_TO/expands/J3_Campbell_heatmap.eps',height = 20)
pheatmap(heat_GBM,na_col = 'gray',drop_levels = T,cluster_rows = FALSE, cluster_cols = FALSE, main = 'J3 - Campbell Cancer Driver Genes',cellheight = 10,cellwidth = 10,color = c('white','red','blue','springgreen4'),breaks = c(0,0.5,1.5,3.5,5))
dev.off()

# J7
heat_GBM = J7[rownames(J7) %in% campbell,]
setEPS()
postscript('/Volumes/BioDiscovery/Lab/Tofilon/P06_Cambridge_Joseph/cellLines_TO/expands/J7_Campbell_heatmap.eps',height = 20)
pheatmap(heat_GBM,na_col = 'gray',drop_levels = T,cluster_rows = FALSE, cluster_cols = FALSE, main = 'J7 - Campbell Cancer Driver Genes',cellheight = 10,cellwidth = 10,color = c('white','red','blue','springgreen4'),breaks = c(0,0.5,1.5,3.5,5))
dev.off()

# J14
heat_GBM = J14[rownames(J14) %in% campbell,]
setEPS()
postscript('/Volumes/BioDiscovery/Lab/Tofilon/P06_Cambridge_Joseph/cellLines_TO/expands/J14_Campbell_heatmap.eps',height = 20)
pheatmap(heat_GBM,na_col = 'gray',drop_levels = T,cluster_rows = FALSE, cluster_cols = FALSE, main = 'J14 - Campbell Cancer Driver Genes',cellheight = 10,cellwidth = 10,color = c('white','red','blue','springgreen4'),breaks = c(0,0.5,1.5,3.5,5))
dev.off()






#### ------ script below not used --------#####

# library(tidyverse)
# library(gridExtra)
# library(grid)
# library(png)
# library(downloader)
# library(grDevices)
# 
# comb2pngs <- function(imgs, bottom_text = NULL){
#   img1 <-  grid::rasterGrob(as.raster(readPNG(imgs[1])),
#                             interpolate = FALSE)
#   img2 <-  grid::rasterGrob(as.raster(readPNG(imgs[2])),
#                             interpolate = FALSE)
#   grid.arrange(img1, img2, nrow = 1, bottom = bottom_text)
# }
# 
# png1_path <- "/Users/valdezkm/Documents/Cambridge/superFreq/myAnalysis/plots/J3/scatters/J3_norm/J3_T1/all.png"
# png2_path <- "/Users/valdezkm/Documents/Cambridge/superFreq/myAnalysis/plots/J3/scatters/J3_norm/J3_T1/allNamed.png"
# 
# comb2pngs(c(png1_path, png2_path))


