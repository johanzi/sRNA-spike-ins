#examine.MFE.distributions.R
#generates violin plots of MFE distributions for endogeneous miRNAs and 252 sets of 65,536 sequences (9 at a time; 28 individual plots) for examination of MFE distributions
#also plots separate graph showing violin plots corresponding to 8 sets of random sequences used in this study

library("extrafont")
library("vioplot")
library("RColorBrewer")

#setwd("/Volumes/compFiles/sRNAs/projects/apps/smallRNA_spikeIns/scripts.for.ms/methods.sRNA.spike.in.design")
setwd("S:/SCRIPTS/sRNA-spike-ins/methods.sRNA.spike.in.design")

#function to plot violins
plot.violins <- function(fileList, name) {
  miRNA.mfes = read.delim("mature.miRNA.seqs.top50percent_mfes", header=FALSE)
  
  nameList = c("ctrl")
  mfeList = c(miRNA.mfes)
  
  for (file in fileList) {
    file.name = strsplit(file,"_")[[1]][1]
    nameList = c(nameList,file.name)
    file.mfes = read.delim(paste("randomOligoSets/folded/mfes/",file,sep=""), header=FALSE)
    mfeList = c(mfeList,file.mfes)
  }
  
  outFile = paste("graphs/",name,".violins.pdf",sep='')
  
  col_v = brewer.pal(10, "Set3")
  
  #pdf(outFile, family='Arial', useDingbats=F) # Error unknown family 'Arial'
  pdf(outFile, family='Helvetica', useDingbats=F)
  vioplot(mfeList[1]$V1,mfeList[2]$V1,mfeList[3]$V1,mfeList[4]$V1,mfeList[5]$V1,mfeList[6]$V1,mfeList[7]$V1,mfeList[8]$V1,mfeList[9]$V1,mfeList[10]$V1, col="white", names=nameList)
  vioplot(mfeList[1]$V1, col=col_v[1], at=1, add=T)
  vioplot(mfeList[2]$V1, col=col_v[2], at=2, add=T)
  vioplot(mfeList[3]$V1, col=col_v[3], at=3, add=T)
  vioplot(mfeList[4]$V1, col=col_v[4], at=4, add=T)
  vioplot(mfeList[5]$V1, col=col_v[5], at=5, add=T)
  vioplot(mfeList[6]$V1, col=col_v[6], at=6, add=T)
  vioplot(mfeList[7]$V1, col=col_v[7], at=7, add=T)
  vioplot(mfeList[8]$V1, col=col_v[8], at=8, add=T)
  vioplot(mfeList[9]$V1, col=col_v[9], at=9, add=T)
  vioplot(mfeList[10]$V1, col=col_v[10], at=10, add=T)
  abline(h=0, lty=2)
  dev.off()
}

fileList = list.files("randomOligoSets/folded/mfes/", pattern = "mfes")
list.1 = fileList[1:9]
list.2 = fileList[10:18]
list.3 = fileList[19:27]
list.4 = fileList[28:36]
list.5 = fileList[37:45]
list.6 = fileList[46:54]
list.7 = fileList[55:63]
list.8 = fileList[64:72]
list.9 = fileList[73:81]
list.10 = fileList[82:90]
list.11 = fileList[91:99]
list.12 = fileList[100:108]
list.13 = fileList[109:117]
list.14 = fileList[118:126]
list.15 = fileList[127:135]
list.16 = fileList[136:144]
list.17 = fileList[145:153]
list.18 = fileList[154:162]
list.19 = fileList[163:171]
list.20 = fileList[172:180]
list.21 = fileList[181:189]
list.22 = fileList[190:198]
list.23 = fileList[199:207]
list.24 = fileList[208:216]
list.25 = fileList[217:225]
list.26 = fileList[226:234]
list.27 = fileList[235:243]
list.28 = fileList[244:251]
#list.28 = fileList[244:252]

#step through list by 9s and plot to examine corresponding MFE distributions
plot.violins(list.1, "Set1")
plot.violins(list.2, "Set2")
plot.violins(list.3, "Set3")
plot.violins(list.4, "Set4")
plot.violins(list.5, "Set5")
plot.violins(list.6, "Set6")
plot.violins(list.7, "Set7")
plot.violins(list.8, "Set8")
plot.violins(list.9, "Set9")
plot.violins(list.10, "Set10")
plot.violins(list.11, "Set11")
plot.violins(list.12, "Set12")
plot.violins(list.13, "Set13")
plot.violins(list.14, "Set14")
plot.violins(list.15, "Set15")
plot.violins(list.16, "Set16")
plot.violins(list.17, "Set17")
plot.violins(list.18, "Set18")
plot.violins(list.19, "Set19")
plot.violins(list.20, "Set20")
plot.violins(list.21, "Set21")
plot.violins(list.22, "Set22")
plot.violins(list.23, "Set23")
plot.violins(list.24, "Set24")
plot.violins(list.25, "Set25")
plot.violins(list.26, "Set26")
plot.violins(list.27, "Set27")
plot.violins(list.28, "Set28")


#plot for random sequence sets used in this study

miRNA.mfes = read.delim("mature.miRNA.seqs.top50percent_mfes", header=FALSE)

nameList = c("ctrl")
mfeList = c(miRNA.mfes)

fileList = c("433_mfes","403_mfes","361_mfes","871_mfes","974_mfes","71_mfes","823_mfes","87_mfes")

for (file in fileList) {
  file.name = strsplit(file,"_")[[1]][1]
  nameList = c(nameList,file.name)
  file.mfes = read.delim(paste("randomOligoSets/folded/mfes/",file,sep=""), header=FALSE)
  mfeList = c(mfeList,file.mfes)
}

outFile = "../supplemental.figure.4/supplemental.figure.4.pdf"

col_v = brewer.pal(9, "Set3")

pdf(outFile, family='Helvetica', useDingbats=F)
vioplot(mfeList[1]$V1,mfeList[2]$V1,mfeList[3]$V1,mfeList[4]$V1,mfeList[5]$V1,mfeList[6]$V1,mfeList[7]$V1,mfeList[8]$V1,mfeList[9]$V1, col="white", names=nameList)
vioplot(mfeList[1]$V1, col=col_v[1], at=1, add=T)
vioplot(mfeList[2]$V1, col=col_v[2], at=2, add=T)
vioplot(mfeList[3]$V1, col=col_v[3], at=3, add=T)
vioplot(mfeList[4]$V1, col=col_v[4], at=4, add=T)
vioplot(mfeList[5]$V1, col=col_v[5], at=5, add=T)
vioplot(mfeList[6]$V1, col=col_v[6], at=6, add=T)
vioplot(mfeList[7]$V1, col=col_v[7], at=7, add=T)
vioplot(mfeList[8]$V1, col=col_v[8], at=8, add=T)
vioplot(mfeList[9]$V1, col=col_v[9], at=9, add=T)
abline(h=0, lty=2)
dev.off()