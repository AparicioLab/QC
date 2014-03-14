# R script to read in the calculated and obtained Tms of amplicons
# Change the sample name between #### in this script to make it run correctly

# install.packages("calibrate")
library(calibrate)

#######################################################################
RUNID="A7DWM"
Project="IntClust"
# if run directly uncomment the sample name
# Command line `Rscript Temp_Corr2.R --no-save --no-restore --args SA499`

# This takes the 4th argument (see str above) which is sample name
# args <- commandArgs(trailingOnly = TRUE)
# SA <- args[4]

# This is where the input files from script are deposited

# This is where the output of this R-script will be
path = paste(paste("C:\\Users\\dyap_000\\Documents\\R",Project,sep="\\"),RUNID,sep="\\")
setwd(path)

#########################################################################

check="freq"
createCounter <- function(value) { function(i) { value <<- value+i} }
count <- createCounter(1)

############# DO NOT CHANGE ANYTHING BELOW THIS LINE #################

# Find the file names with the same name and pattern in the working dir
pat=paste("*.*","allfreq.csv",sep="-")
file_names = list.files(pattern = pat)

# Reads the file name into Q (in a list)
QC_list = lapply(file_names, read.csv, header = TRUE)
# Change the number of samples

# check with "list(QC_list)"
conditions <- length(QC_list)

# assume all have same no of samples
sample1 <- ncol(QC_list[[1]])
sample2 <- ncol(QC_list[[2]])

if (sample1 == sample2)	samples <- sample1

len1<-nrow(QC_list[[1]])
len2<-nrow(QC_list[[2]])

if (len1 > len2)	len <- len1
if (len1 < len2)	len <- len2
if (len1 == len2)	len <- len1

d.frame <- data.frame(	     ID = rep("", len),
 			     FFPE = rep (0, len),
 			     Frozen = rep (0, len),
			     stringsAsFactors = FALSE)
# Col1 is the ID
for (rj in 2:samples ) 	{

sum1<-"NULL"
assign(paste("Sample", rj, sep=""), sum1)

FFPE<-QC_list[[1]][c(1,rj)]
FFPEa <- na.omit(FFPE[order(FFPE$ID), ]) 
FFPEdat <- FFPEa[which(!duplicated(FFPEa$ID)),]

Frozen<-QC_list[[2]][c(1,rj)]
Frozena <- na.omit(Frozen[order(Frozen$ID), ]) 
Frozendat <- Frozena[which(!duplicated(Frozena$ID)),]

sum1 <- merge(FFPEdat, Frozendat, by="ID", all=TRUE)

assign(paste("Sample", rj, sep=""), sum1)


        name <- paste("Frozen-FFPE-Sample", rj, sep="")
        fname <- paste(path,name,sep="\\")
        

        # Output each as a separate PDF
        pdf (file=paste(fname, "pdf", sep="."))

        # Formats the plots
        title <- paste("Correlation of Methods for", name, sep=" ")
        plot(sum1, title=title)
        dev.off()

	csvfile<- paste(fname,"csv",sep=".")
	write.table(sum1,file=csvfile,sep=",",row.names=FALSE,col.names=TRUE)

 

        print("end")

}


