#examine.MFE.distributions.R
#generates violin plots of MFE distributions for endogeneous miRNAs and 252 sets of 65,536 sequences (9 at a time; 28 individual plots) for examination of MFE distributions
#also plots separate graph showing violin plots corresponding to 8 sets of random sequences used in this study

library("extrafont")
library("RColorBrewer")
library("vioplot")



# Variables 

# Output directory
setwd("S:/SCRIPTS/sRNA_spike_in_design/maize")

# Name of the MFE file for miRNAs
miRNA_file <- "mature.miRNA.seqs_mfes"

# Location of MFE files
location_mfes <- "randomOligoSets/folded/mfes/"


# Output directory for the violin graphs
location_graphs <-  "graphs/"



#function to plot violins

plot.violins <- function(location_graphs, location_mfes, fileList, name) {
  
  # Read content of MFE for miRNAs
  miRNA.mfes = read.delim(miRNA_file, header=FALSE)
  
  # Create a name and MFE list, containing first name and data for miRNAs
  nameList = c("miRNAs")
  mfeList = c(miRNA.mfes)

  # Loop over the MFE files to retrieve names and append data to mfeList and nameList
  for (file in fileList) {
    file.name = strsplit(file,"_")[[1]][1]
    nameList = c(nameList, file.name)
    file.mfes = read.delim(paste(location_mfes, file, sep=""), header=FALSE)
    mfeList = c(mfeList, file.mfes)
    names
  }
  
  # Generate output PDF file
  outFile = paste(location_graphs, name, "_violins.pdf", sep='')
  
  # Define colors (as many as there are files in the list
  col_v = brewer.pal(11, "Set3")
  col_v_used <- col_v[1:length(nameList)]
  
  # Open PDF file to write
  pdf(outFile, family='Helvetica', useDingbats=F)
  
  # Plot
  vioplot(mfeList, names=as.vector(nameList), cex.axis = 0.6, col=col_v_used)
 
  # Add a dashed line at y=0
  abline(h=0, lty=2)
  dev.off()
  
}


# Create a list of list of MFE files

create_list_files <- function(location_mfes){

  # Get a list of all files containing suffix "mfes"
  fileList = list.files(location_mfes, pattern = "mfes")
  
  # Check number of files
  len_list <- length(fileList)
  
  # Create a list of each file and pool them 10 by 10
  # Initialize the first element of the list
  y=1
  # Create empty receiving list
  list_list_files <- list()
  
  for(i in seq(0,len_list,10)){
    start=i+1
    end=i+10
    list_list_files[[y]] <- fileList[start:end]
    y=y+1
  }
  
  # Remove any NAs from the lists (usually in the last list if not a multiple of 10)
  list_list_files_filtered <- lapply(list_list_files, function(x) x[!is.na(x)])
  
  # Return curated list
  return(list_list_files_filtered)
  
}

# Create of list of 10 by 10 mfe files
list_list_files <- create_list_files(location_mfes)


#step through list by 10 and plot to examine corresponding MFE distributions
for(i in (1:length(list_list_files))){
  name_plot = paste("Set",i ,sep="")
  plot.violins(location_graphs, location_mfes, list_list_files[[i]], name_plot)
}


# Test
plot.violins(location_graphs, location_mfes, list_list_files[[26]], "test")





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