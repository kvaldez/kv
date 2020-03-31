######### Gene level QC and DEG Prep - TMM Normalization #########
#### Start prepare data for DEG analysis ####
library(reshape2)
library(ggplot2)
library(edgeR)
library(plotly)
library(stringr)
library(purrr)
library(tidyverse)




#### Human Filtered TMM normalization ####
load('/Volumes/BioDiscovery/Lab/Tofilon/P16_OB_RH/Counted/fc_humanFiltered.RData')
#
geneCounts = data.frame(fc_humanFiltered$counts)
colnames(geneCounts) = gsub('\\.bam','',colnames(geneCounts))
#
sampleinfo = data.frame(sampleNames = colnames(geneCounts))
sampleinfo$condition = c(rep('OB',4),rep('RH',4))
sampleinfo$label = gsub('\\..*','',sampleinfo$sampleNames)  
sampleinfo$sampleNames == colnames(geneCounts)
sampleinfo$mouse = str_sub(sampleinfo$label,start = -1)  
sampleinfo$mouse = paste0('M',sampleinfo$mouse)
#
geneCounts = DGEList(counts=geneCounts)
keep <- rowSums(cpm(geneCounts)>0.5) >= 2 #                                         #filter lowly expressed genes
geneCounts <- geneCounts[keep,,keep.lib.sizes = FALSE]
geneCounts = calcNormFactors(geneCounts,method='TMM')     #                         #normalize with TMM
#
geneCountsCPM = cpm(geneCounts,log = TRUE, normalized.lib.sizes = TRUE, prior.count = 0.5)#     Norm lib sizes TRUE
# write.csv(geneCountsCPM,'/Volumes/BioDiscovery/Lab/Tofilon/P16_OB_RH/Counted/normalized_log2_CPM_TMM.csv')
#### End prepare data for DEG analysis ####

#
#####  Post normalization QC  #####
df.m = melt(geneCountsCPM)
#
# Boxplot
par(mar=c(par("mar")[1]+5,par("mar")[-1]))
boxplot(value~Var2,las=2,data=df.m,main="Filtered Counts After TMM Normalization",
        ylab="Log2 CPM",col=as.numeric(as.factor(sampleinfo$condition))*5,
        names=sampleinfo$label)

# Histogram
beforehist <- ggplotly(ggplot(df.m) + geom_density(aes(x = value,colour = Var2)) + labs(x = NULL) + theme(legend.position='right') + scale_x_log10())
beforehist

#PCA plotly
before.edf=geneCountsCPM
before.tedf= t(before.edf)
before.pca=prcomp(before.tedf,scale.=T)
Phenotype=sampleinfo$condition
cell_rep=sampleinfo$label
before.pc1 = round(before.pca$sdev[1]^2/sum(before.pca$sdev^2)*100,2)
before.pc2 = round(before.pca$sdev[2]^2/sum(before.pca$sdev^2)*100,2)
before.pc3 = round(before.pca$sdev[3]^2/sum(before.pca$sdev^2)*100,2)
pcafactor = as.factor(sampleinfo$condition)
library(RColorBrewer)
col <- brewer.pal(nlevels(pcafactor), "Paired")
p <- plot_ly(as.data.frame(before.pca$x[,1:3]), x = ~PC1, y = ~PC2, z = ~PC3, color = pcafactor, colors = col, hoverinfo="text",
             text = ~sampleinfo$label,
             hovertext = ~sampleinfo$label, textposition = 'middle right',
             textfont = list(color = '#000000', size = 8))  %>%
  add_markers() %>%
  layout(title = "After TMM Normalization PCA plot",
         scene = list(xaxis = list(title = paste0("PC1 (",before.pc1,"%)")),
                      yaxis = list(title = paste0("PC2 (",before.pc2,"%)")),
                      zaxis = list(title = paste0("PC3 (",before.pc3,"%)"))))
p

# PCA rgl
library(rgl)
open3d()
pal=colorRamps::primary.colors(19)
#pal=colorRamps::primary.colors(2)
plot3d(before.pca$x[,1:3],
       xlab = paste0("PC1 (",before.pc1,"%)"),
       ylab = paste0("PC2 (",before.pc2,"%)"),
       zlab = paste0("PC3 (",before.pc3,"%)"),
       col=pal[as.numeric(as.factor(sampleinfo$condition))], type='s',size=2)
group.v=as.vector(sampleinfo$label)
text3d(before.pca$x, before.pca$y, before.pca$z, group.v, cex=0.6, adj=1.5)
dev.off()

# Similarity Heatmap
library(amap)
library(pheatmap)
mat=as.data.frame(as.matrix(Dist(t(geneCountsCPM),method = 'pearson',diag = TRUE)))
mat = 1 - mat
pheatmap(mat,labels_col = as.character(sampleinfo$label), labels_row = as.character(sampleinfo$label))



# Anova function
library(gtools)
set.seed(123)
posthoc_withMean <- function(dat,Categories) {
  pVal.posthoc = matrix(NA, nrow = nrow(dat), ncol = (length(c(combinations(length(unique(Categories)),2)))+length(unique(Categories))))
  pval = c()
  for(i in 1: nrow(dat)){
    fm1 = aov(dat[i,]~ Categories) #            
    tmp = TukeyHSD(fm1, "Categories", ordered = F)$Categories
    tmp1 = matrix(NA, ncol = ncol(pVal.posthoc))
    lab = c(paste(rownames(tmp),colnames(tmp)[1],sep="_"),paste(rownames(tmp),colnames(tmp)[4],sep="_"),paste0(unique(Categories),'_mean'))
    means = c()
    for (j in 1:length(unique(Categories))) {
      tempMean = mean(dat[i,Categories %in% unique(Categories)[j]])
      means = c(means,tempMean)
    }
    tmp1 = c(tmp[,1],tmp[,4],means)
    pVal.posthoc[i,] = tmp1
    pval = c(pval,summary(fm1)[[1]][4:5][[2]][1])
  }
  colnames(pVal.posthoc) = lab
  rownames(pVal.posthoc) = rownames(dat)
  pVal.posthoc = data.frame("pval" = round(pval,3),pVal.posthoc)
  pVal.posthoc = round(pVal.posthoc,3)
  #pVal.posthoc = pVal.posthoc[,grep("pval|p.adj",colnames(pVal.posthoc))]#to get only pvals
  return(pVal.posthoc)
}



#### ANOVA Post Hoc DEG - OB vs RH ####
dat = geneCountsCPM
colnames(dat) == sampleinfo$sampleNames#      Check to make sure names match pheno
colnames(dat) = sampleinfo$label
Categories = sampleinfo$condition
#
res1 = posthoc_withMean(dat,Categories)          # only run again if nec, takes forever
tmp = res1
tmp$Symbol = gsub('\\..*','',rownames(tmp))
tmp$Ensembl = rownames(tmp)
annot_ensml = read.delim('/Users/valdezkm/Documents/annotation_Ensembl.txt',sep = '\t')
tmp = merge(tmp,annot_ensml,by.x='Symbol',by.y='Gene.stable.ID',all.x=T)
tmp = tmp[order(tmp$pval),]
tmp = tmp[,c(8,3:7,9,11)]
sig = tmp[tmp$RH.OB_p.adj<0.05,]
#
sig_up = sig[sig$RH.OB_diff>=0,]
sig_down = sig[sig$RH.OB_diff<0,]
#
write.table(tmp,'/Volumes/BioDiscovery/Lab/Tofilon/P16_OB_RH/ANOVA/DEG.txt',row.names = F,sep='\t')





######## Heatmap of Venn Genes ##########
library(pheatmap)
heat = geneCountsCPM
# batch = ifelse(grepl('2|6',sampleinfo$label),'batch_2','batch_1')
#
# modcombat = model.matrix(~as.factor(condition), data=sampleinfo)
# heat = ComBat(dat = heat, batch = batch, mod = modcombat, par.prior = TRUE, prior.plots = TRUE)
heat = heat[rownames(heat) %in% sig$Ensembl,]
heat = t(scale(t(heat)))
#
pheatmap(heat,drop_levels = T, labels_col = as.character(sampleinfo$label),fontsize_row = 1,clustering_distance_cols = 'correlation',clustering_distance_rows = 'correlation')







#### Individual Mouse DEG #####
sampleinfo = data.frame(sampleNames = colnames(geneCounts))
sampleinfo$condition = c(rep('OB',4),rep('RH',4))
sampleinfo$label = gsub('\\..*','',sampleinfo$sampleNames)  
sampleinfo$sampleNames == colnames(geneCounts)
sampleinfo$mouse = str_sub(sampleinfo$label,start = -1)  
sampleinfo$mouse = paste0('M',sampleinfo$mouse)
#
dat = geneCountsCPM
colnames(dat) == sampleinfo$sampleNames#      Check to make sure names match pheno
colnames(dat) = sampleinfo$label
#
# subset data based on mouse
M2 = which(sampleinfo$mouse=='M2')
M4 = which(sampleinfo$mouse=='M4')
M6 = which(sampleinfo$mouse=='M6')
M7 = which(sampleinfo$mouse=='M7')
#
allDat = list(M2=dat[,M2], M4=dat[,M4], M6=dat[,M6], M7=dat[,M7])

# simple subtraction of gene log2(cpm) between RH-OB for each mouse
for (i in 1:length(allDat)) {
  allDat[[i]] = as.data.frame(allDat[[i]])
  allDat[[i]]$diff = allDat[[i]][,2] - allDat[[i]][,1]
  allDat[[i]]$Symbol = gsub('\\..*','',rownames(allDat[[i]]))
  allDat[[i]]$Ensembl = rownames(allDat[[i]])
  annot_ensml = read.delim('/Users/valdezkm/Documents/annotation_Ensembl.txt',sep = '\t')
  allDat[[i]] = merge(allDat[[i]],annot_ensml,by.x='Symbol',by.y='Gene.stable.ID',all.x=T)
  allDat[[i]] = allDat[[i]][,c(6,2:5,7,9)]
  allDat[[i]] = allDat[[i]][order(abs(allDat[[i]]$diff),decreasing = TRUE),]
  #
  # write.table(allDat[[i]],paste0('/Volumes/BioDiscovery/Lab/Tofilon/P16_OB_RH/Mouse/',names(allDat[i]),'.txt'),row.names = F,sep = '\t')
}

# Venn prep: |logFC| > 2
mouseUp = allDat
mouseDown = allDat
for (i in 1:length(mouseUp)) {
  mouseUp[[i]] = mouseUp[[i]][mouseUp[[i]]$diff>=2,]
  mouseDown[[i]] = mouseDown[[i]][mouseDown[[i]]$diff<=(-2),]
}

VennFour_expr = function(first,nameFirst,second,nameSecond,third,nameThird,fourth,nameFourth) {
  library(gdata)
  library(gplots)
  library(stringr)
  
  first = as.data.frame(first)
  first = as.list(as.data.frame(first$Ensembl))
  
  second = as.data.frame(second)
  second = as.list(as.data.frame(second$Ensembl))
  
  third = as.data.frame(third)
  third = as.list(as.data.frame(third$Ensembl))
  
  fourth = as.data.frame(fourth)
  fourth = as.list(as.data.frame(fourth$Ensembl))
  
  venn_list = cbind(first, second, third, fourth)
  venn_list <- lapply(as.list(venn_list), function(x) x[x != ""])
  names(venn_list) = c(nameFirst,nameSecond,nameThird,nameFourth)
  inEach = venn(venn_list)
  inters = attr(inEach, "intersections")
  
  inters
}

vUp = VennFour_expr(mouseUp$M2,'M2',mouseUp$M4,'M4',mouseUp$M6,'M6',mouseUp$M7,'M7')
vDown = VennFour_expr(mouseDown$M2,'M2',mouseDown$M4,'M4',mouseDown$M6,'M6',mouseDown$M7,'M7')



#### Create Venn TABLE ####
subUp = mouseUp
for (i in 1:length(mouseUp)) {
   subUp[[i]] = subUp[[i]][,c(1,5,6,7,2,3,4)]
   colnames(subUp[[i]])[colnames(subUp[[i]])=='diff'] = paste0(names(subUp[i]),"_diff")
}
subDown = mouseDown
for (i in 1:length(mouseDown)) {
  subDown[[i]] = subDown[[i]][,c(1,5,6,7,2,3,4)]
  colnames(subDown[[i]])[colnames(subDown[[i]])=='diff'] = paste0(names(subDown[i]),"_diff")
}
tableUp = subUp %>% reduce(full_join,by=colnames(subUp[[1]])[1:4])
tableDown = subDown %>% reduce(full_join,by=colnames(subDown[[1]])[1:4])
tableUp$direction = 'Up'
tableDown$direction = 'Down'

vTable = rbind(tableUp,tableDown)
write.table(vTable,'/Volumes/BioDiscovery/Lab/Tofilon/P16_OB_RH/Mouse/Venn_Table.txt',sep='\t',row.names = F, na = '')








