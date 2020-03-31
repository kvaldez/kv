library(Hmisc)
library(rlist)
library(tidyverse)
library(purrr)
library(amap)
library(pheatmap)
library(plotly)
library(ggplot2)
library(reshape2)


#### TUMOR ONLY: A25: Read in mafs, filter -> Vaf heatmap ####
temp = list.files(path = '/Volumes/BioDiscovery/Lab/Tofilon/P07_A25_J14/A25/mutect2',pattern = '*.maf',full.names = TRUE)
allMafs = lapply(temp, read.delim, stringsAsFactors = F)
names(allMafs) = gsub('S._|_in.*','',basename(temp))

# add vaf
for (i in 1:length(allMafs)) {
  allMafs[[i]]$vaf = allMafs[[i]]$t_alt_count / allMafs[[i]]$t_depth
  allMafs[[i]] = allMafs[[i]][allMafs[[i]]$vaf > 0.05,]
}

# filter out non-protein-altering variants
variantsToKeep = c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site","Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del","In_Frame_Ins", "Missense_Mutation")
filtMafs = list()
for (i in 1:length(allMafs)) {
  filtMafs[[i]] = allMafs[[i]][allMafs[[i]]$Variant_Classification %in% variantsToKeep,]
  names(filtMafs)[i] = names(allMafs)[i]
}

# # write out filtered variants
# colsToKeep = c(1:19,35:76,133)
# temp = filtMafs
# for (i in 1:length(filtMafs)) {
#   temp[[i]] = temp[[i]][,colsToKeep]
#   write.csv(temp[[i]],paste0('/Volumes/BioDiscovery/Lab/Tofilon/P07_A25_J14/A25/mutect2/filtered_variants/',names(temp)[i],'_filteredVariants.csv'),row.names = F)
# }



#### Pearsons correlation coefficients - VAFs - filtered data ####
justVAFs = c()
for (i in 1:length(filtMafs)) {
  justVAFs[[i]] = subset(filtMafs[[i]],select=c(Hugo_Symbol,Chromosome,Start_Position,vaf))
  colnames(justVAFs[[i]])[4] = names(filtMafs[i])
}
df_justVAFs = justVAFs %>% reduce(full_join, by = c('Hugo_Symbol','Chromosome','Start_Position'))
mat_justVAFs = as.matrix(df_justVAFs[,4:ncol(df_justVAFs)])

#
# heatmap - pearson correlation using rcorr
cors = rcorr(mat_justVAFs, type = 'pearson')
cors_r = cors$r
pheatmap(cors_r,display_numbers = T,fontsize = 14)



#### Subset GBM driver genes ####
JoeGenes = read.delim('/Volumes/BioDiscovery/Lab/Tofilon/Joseph_GBM_Drivers.txt')
JoeGenes = as.character(unique(JoeGenes$GBM.Drive.Genes))
# subset driver genes
gbmMafs = list()
for (i in 1:length(filtMafs)) {
  gbmMafs[[i]] = filtMafs[[i]][filtMafs[[i]]$Hugo_Symbol %in% JoeGenes,]
  names(gbmMafs)[i] = names(filtMafs)[i]
}
# heatmap of GBM variants and vafs
GBMvar = list()
for (i in 1:length(gbmMafs)) {
  if (nrow(gbmMafs[[i]]) > 0) {
    GBMvar[[i]] = subset(gbmMafs[[i]],select=c(Hugo_Symbol,HGVSc,vaf))
    GBMvar[[i]]$geneVar = paste(GBMvar[[i]]$Hugo_Symbol,GBMvar[[i]]$HGVSc)
    GBMvar[[i]]$Sample = names(gbmMafs)[i]
  } else {
    GBMvar[[i]] = data.frame(matrix(ncol=3))
    colnames(GBMvar[[i]]) = c('geneVar','vaf')
    GBMvar[[i]]$geneVar = as.character(GBMvar[[i]]$geneVar)
    GBMvar[[i]]$vaf = as.numeric(GBMvar[[i]]$vaf)
    GBMvar[[i]]$Sample = names(gbmMafs)[i]
  }
  GBMvar[[i]] = subset(GBMvar[[i]],select=c(geneVar,vaf,Sample))
  names(GBMvar)[i] = names(gbmMafs)[i]
}

df_GBMvar = GBMvar %>% reduce(full_join,by=c('geneVar','vaf','Sample'))
df_GBMvar = dcast(df_GBMvar,geneVar~Sample,value.var = 'vaf')
df_GBMvar = df_GBMvar[!is.na(df_GBMvar$geneVar),]
rownames(df_GBMvar) = df_GBMvar$geneVar
df_GBMvar = df_GBMvar[,-1]
textVar = df_GBMvar
textVar = round(textVar,digits = 2)
textVar[is.na(textVar)] = ''
# Edit long variant name
rownames(df_GBMvar)[10] = gsub('ins.*','\\3',rownames(df_GBMvar)[10])

pheatmap(df_GBMvar,treeheight_row = 0,treeheight_col = 0,display_numbers = textVar, fontsize = 12, cluster_cols = F,cluster_rows = FALSE, na_col = 'white')



#### Copy number correlations ####
# read in CNVs
temp = list.files(path = '/Volumes/BioDiscovery/Lab/Tofilon/P07_A25_J14/A25/cnvkit',pattern = '*_calls.cns',full.names = TRUE)
allCNVs = lapply(temp, read.delim, stringsAsFactors = F)
names(allCNVs) = gsub('S._|_in.*','',basename(temp))

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
      if ((snv$Chromosome[i] == cnv$chromosome[j]) & (snv$Start_Position[i] >= cnv$start[j]) & (snv$End_Position[i] <= cnv$end[j])) {
        snv$CNt[i] = cnv$cn[j]
        snv$A[i] = cnv$cn1[j]
        snv$B[i] = cnv$cn2[j]
      }
    }
  }
  return(snv)
}

SNV_CNV = c()
for (k in 1:length(allCNVs)) {
  SNV_CNV[[k]] = assignCNVtoSNV(filtMafs[k],allCNVs[k])
}
names(SNV_CNV) = names(filtMafs)

justCNs = c()
for (i in 1:length(SNV_CNV)) {
  justCNs[[i]] = subset(SNV_CNV[[i]],select=c(Hugo_Symbol,Chromosome,Start_Position,CNt))
  colnames(justCNs[[i]])[4] = names(SNV_CNV[i])
}

df_justCNs = justCNs %>% reduce(full_join, by = c('Hugo_Symbol','Chromosome','Start_Position'))
mat_justCNs = as.matrix(df_justCNs[,4:ncol(df_justCNs)])

cors = rcorr(mat_justCNs, type = 'pearson')
cors_r = cors$r
pheatmap(cors_r,display_numbers = T,fontsize = 14)


#### Create filtered mafs for expands input (too many variants for expands to work properly) ####
expandsMafs = filtMafs
names(expandsMafs) = gsub('.maf','',basename(temp))  
for (i in 1:length(expandsMafs)) {
  write.table(expandsMafs[[i]],paste0('/Volumes/BioDiscovery/Lab/Tofilon/P07_A25_J14/A25/expands_filteredMafs/',names(expandsMafs[i]),'_filtered.maf'),quote = F,sep = '\t',row.names = F)
}


#### Expands - read in number of SPs per sample ####
temp = list.files(path = '/Volumes/BioDiscovery/Lab/Tofilon/P07_A25_J14/A25/expands_output',pattern = '*dm.csv',full.names = TRUE)
expands = lapply(temp, read.delim, stringsAsFactors = F, sep=',')
names(expands) = gsub('S._|_in.*','',basename(temp))
#
#
numSPs = data.frame(Sample = c(names(expands)))
for (i in 1:length(expands)) {
  numSPs$numSPs[[i]] = nrow(unique(expands[[i]][12]))
}
numSPs = numSPs[order(numSPs$Sample),]
# write.csv(numSPs,'/Volumes/BioDiscovery/Lab/Tofilon/P07_A25_J14/A25/expands_output/numSPs.csv',row.names = F)


#### Create expands subclone plots ####
####merge expands output with mafs ####
# GBM driver genes
# JoeGenes = read.delim('/Volumes/BioDiscovery/Lab/Tofilon/Joseph_GBM_Drivers.txt')
# JoeGenes = as.character(unique(JoeGenes$GBM.Drive.Genes))
#
# Campbell Cancer Driver Genes
campbell = read.delim('/Volumes/BioDiscovery/Lab/Tofilon/P06_Cambridge_Joseph/Campbell_Cancer_Genes.csv',sep=',',header=F)
campbell = unlist(campbell)
campbell = campbell[campbell!='']
campbell = droplevels(campbell)

#
# sort mafs and expands lists so they match
filtMafs = filtMafs[order(names(filtMafs))]
expands = expands[order(names(expands))]
names(filtMafs) == names(expands)
#
# annotate expands..... if all variants run through exands, keep all filtered mafs, omit non-protein-altering and low vafs with all.x=TRUE
expands_annot = list()
for (i in 1:length(expands)) {
  filtMafs[[i]]$Chromosome = gsub('chr','',filtMafs[[i]]$Chromosome)
  expands_annot[[i]] = merge(filtMafs[[i]],expands[[i]],by.x=c('Chromosome','Start_Position'),by.y=c('chr','startpos'),all.x=TRUE)
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
expands_annot = expands_annot %>% reduce(full_join)
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
expands_annot = heat(expands_annot)
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
heat_GBM = expands_annot[rownames(expands_annot) %in% campbell,]
setEPS()
postscript('/Volumes/BioDiscovery/Lab/Tofilon/P07_A25_J14/A25/expands_output/Campbell_heatmap.eps')
pheatmap(heat_GBM,na_col = 'gray',drop_levels = T,cluster_rows = FALSE, cluster_cols = FALSE, main = 'Campbell Cancer Driver Genes',cellheight = 10,cellwidth = 10,color = c('white','red','blue','springgreen4'),breaks = c(0,0.5,1.5,3.5,5))
dev.off()






#### YAPSA Cosmic mutational signatures ####
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

# load correlation corrections for later - using WES instead of WGS
# instructions for creating correls found in /Volumes/BioDiscovery/Lab/Tofilon/P06_Cambridge_Joseph/Cambridge.R
load('/Volumes/BioDiscovery/Lab/Tofilon/YAPSA/correls.RData')


####  prepare samples
yap_name = 'A25'
yap = filtMafs
for (i in 1:length(yap)) {
  yap[[i]]$Sample = names(yap[i])
  yap[[i]] = yap[[i]][yap[[i]]$Variant_Type=='SNP',]
  yap[[i]] = yap[[i]][which(colnames(yap[[i]])%in%c('Chromosome','Start_Position','Reference_Allele','Tumor_Seq_Allele2','Sample','chrStart'))]
}
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
  in_subgroups_df = COSMIC_subgroups_df)   

# And relative exposures:
exposures_barplot(
  in_exposures_df = CosmicValid_cutoffSpec_LCDlist$norm_exposures,
  in_signatures_ind_df = CosmicValid_cutoffSpec_LCDlist$out_sig_ind_df,
  in_subgroups_df = COSMIC_subgroups_df)

# cluster samples based on their signature exposures
complex_heatmap_exposures(CosmicValid_cutoffSpec_LCDlist$norm_exposures,
                          COSMIC_subgroups_df,
                          CosmicValid_cutoffSpec_LCDlist$out_sig_ind_df,
                          in_data_type="norm exposures",
                          in_subgroup_colour_column="col",
                          in_method="manhattan",
                          in_subgroup_column="subgroup")






#### ABSOLUTE with filtered data (only protein-coding) - prepare files ####
# upload segmetrics files
temp = list.files(path = '/Volumes/BioDiscovery/Lab/Tofilon/P07_A25_J14/A25/cnvkit',pattern = '*.segmetrics',full.names = TRUE)
segs = lapply(temp, read.delim, stringsAsFactors = F)
names(segs) = gsub('S._|_in.*','',basename(temp))


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
  write.table(segs[[i]],paste0('/Volumes/BioDiscovery/Lab/Tofilon/P07_A25_J14/A25/ABSOLUTE/CN/',names(segs)[i],'.abs.cn.txt'),sep='\t',quote=F,row.names=F)
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
  write.table(abs_mafs[[i]],paste0('/Volumes/BioDiscovery/Lab/Tofilon/P07_A25_J14/A25/ABSOLUTE/MAF/',names(abs_mafs)[i],'.abs.maf.txt'),sep='\t',quote=F,row.names=F,na='')
}

##### Read in ABSOLUTE output ######
temp = list.files(path = '/Volumes/BioDiscovery/Lab/Tofilon/P07_A25_J14/A25/ABSOLUTE/Final_Review',pattern = '_ABS_MAF.txt',full.names = TRUE)
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
abs = abs %>% purrr::reduce(full_join, by = c('genePos','Sample','Pop','Format','Ploidy','cancer_cell_frac'))
abs = dcast(abs,Sample+Pop+Ploidy+Format~genePos,value.var = 'cancer_cell_frac')
abs = as.data.frame(abs)
#
abs = stamppConvert(abs,type='r')

nei = data.frame(stamppNeisD(abs,FALSE))
colnames(nei) = rownames(nei)
pheatmap(nei,display_numbers = T,fontsize = 14)
















