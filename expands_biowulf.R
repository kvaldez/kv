library(expands)
library(phylobase)

# args[1] snv maf
# args[2] cn calls
# args[3] output file prefix
args <- commandArgs(trailingOnly=TRUE)
args_1 = read.delim(args[1],stringsAsFactors = F)
args_2 = read.delim(args[2],stringsAsFactors = F)



create_snv = function(maf) {
  data(snv)
  data(cbs)
  maf$vaf = maf$t_alt_count / maf$t_depth
  maf$PN_B = 0
  maf = maf[maf$Variant_Type=='SNP',]
  maf = subset(maf,select=c(Chromosome,Start_Position,End_Position,Reference_Allele,Tumor_Seq_Allele2,vaf,PN_B))
  colnames(maf) = colnames(snv)
  maf[,4:5][maf[,4:5] == 'A'] = 65
  maf[,4:5][maf[,4:5] == 'C'] = 67
  maf[,4:5][maf[,4:5] == 'G'] = 71
  maf[,4:5][maf[,4:5] == 'T'] = 84
  maf = maf[maf$chr!='X',]
  maf = maf[maf$chr!='Y',]
  maf = sapply(maf,as.numeric)
  maf
}
create_cbs = function(cn) {
  cn = subset(cn,select=c(chromosome,start,end,cn))
  colnames(cn) = colnames(cbs)
  cn = cn[cn$chr!='X',]
  cn = cn[cn$chr!='Y',]
  cn = sapply(cn,as.numeric)
}

snv = create_snv(args_1)
cbs = create_cbs(args_2)

expanded = runExPANdS(snv,cbs)
write.csv(expanded$finalSPs,paste0(args[3],'.SPs.csv'),row.names = F)
write.csv(expanded$dm,paste0(args[3],'.dm.csv'),row.names = F)
save(expanded,file=paste0(args[3],'.all.RData'))
write.csv(expanded$densities,paste0(args[3],'.densities.csv'),row.names = F)
write.csv(expanded$sp_cbs,paste0(args[3],'.sp_cbs.csv'),row.names = F)


#### Merge DM list with mutations ####
annotated = merge(expanded$dm,args_1,by.y=c('Chromosome','Start_Position'),by.x=c('chr','startpos'),all.x=TRUE)
write.csv(annotated,paste0(args[3],'.annotatedSPs.csv'),row.names = F)
