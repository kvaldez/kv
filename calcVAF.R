# read in commented part of vcf (start with at least one #)
openFile = file('/Users/valdezkm/Documents/LP/Val_ClinOmics_Compass_026_T1D_E.strelka.snvs.raw.vcf',open = 'r')
comm = readLines(openFile)
comm = comm[grepl('^#',comm)]
close(openFile)

# read in uncommented vcf table, convert relevant columns to char
dat = read.delim('/Users/valdezkm/Documents/LP/Val_ClinOmics_Compass_026_T1D_E.strelka.snvs.raw.vcf',comment.char = '#', header = FALSE)
dat$V11 = as.character(dat$V11)
dat$V5 = as.character(dat$V5)
dat$V8 = as.character(dat$V8)
dat$VAF = NA

for (i in 1:nrow(dat)) {
  # parse out the number of tumor ACGTs 
  As = as.integer(unlist(strsplit(unlist(strsplit(dat[i,]$V11[1],'\\:'))[5],',')))[1]
  Cs = as.integer(unlist(strsplit(unlist(strsplit(dat[i,]$V11[1],'\\:'))[6],',')))[1]
  Gs = as.integer(unlist(strsplit(unlist(strsplit(dat[i,]$V11[1],'\\:'))[7],',')))[1]
  Ts = as.integer(unlist(strsplit(unlist(strsplit(dat[i,]$V11[1],'\\:'))[8],',')))[1]
  
  # determine depth dependent upon REF variant
  if (dat[i,]$V4=='A') {
    DEPTH = As
  } else if (dat[i,]$V4=='C') {
    DEPTH = Cs
  } else if (dat[i,]$V4=='G') {
    DEPTH = Gs
  } else if (dat[i,]$V4=='T') {
    DEPTH = Ts
  }
  
  # populate VAF column depent upon ALT variant
  dat[i,]$VAF[dat[i,]$V5=='A'] = As / (As + DEPTH)
  dat[i,]$VAF[dat[i,]$V5=='C'] = Cs / (Cs + DEPTH)
  dat[i,]$VAF[dat[i,]$V5=='G'] = Gs / (Gs + DEPTH)
  dat[i,]$VAF[dat[i,]$V5=='T'] = Ts / (Ts + DEPTH)
  
  # recreate info column to include vaf
  dat[i,]$V8 = paste0(dat[i,]$V8,';VAF=',round(dat[i,]$VAF,2))
}
# delete vaf column
dat = subset(dat,select=-VAF)

# write out meta data, append table
write.table(comm,'~/Documents/LP/vaf.vcf',sep='\t',row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(dat,'~/Documents/LP/vaf.vcf',sep='\t',append = TRUE,row.names = FALSE,quote = FALSE,col.names = FALSE)
