# R script to read in the calculated and obtained Tms of amplicons
# Change the sample name between #### in this script to make it run correctly

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

############# DO NOT CHANGE ANYTHING BELOW THIS LINE #################

# Find the file names with the same name and pattern in the working dir
pat=paste("*.*","allfreq.csv",sep="-")
file_names = list.files(pattern = pat)

# Reads the file name into Q (in a list)
QC_list = lapply(file_names, read.csv, header = TRUE)
# Change the number of samples

# check with "list(QC_list)"
samples <- length(QC_list)

createCounter <- function(value) { function(i) { value <<- value+i} }
count <- createCounter(1)


##################################################
# THIS IS THE CURRENT  MODULE THAT WORKS !!!!!   #
##################################################

sumdf <- data.frame(	SNV_ID = rep("", samples),
			Variants = rep(0, samples),
			Rows = rep("", samples),
			stringsAsFactors = FALSE)

for (rj in seq(samples)) 	{
			sid <- colnames(QC_list[[rj]])
			varn <- nrow(QC_list[[rj]])

sumdf$SNV_ID[rj] <- sid
sumdf$Variants[rj] <- varn


# For each of the list of samples
len <- nrow(QC_list[[rj]])

	if ( len > 0 ){


d.frame <- data.frame(	     ID = rep("", len),
 			     framename = rep (0, len),
			     stringsAsFactors = FALSE)

if (check == "freq" ) names(d.frame)[2] <- "Varalfreq"
if (check == "depth" ) names(d.frame)[2] <- "Seqdepth"

for (ri in seq(len) ) 	{

d.frame$ID[ri] <- as.character(QC_list[[rj]][ri])
if (check == "calls")	d.frame$Calls[ri] <- geno(vcf_list[[rj]][ri])$GT
if (check == "freq")	d.frame$Varalfreq[ri] <- geno(vcf_list[[rj]][ri])$VF
if (check == "depth")	d.frame$Seqdepth[ri] <- geno(vcf_list[[rj]][ri])$DP

			} 

# unique freq val or seq depth for summary plot
# must be the same length as names(d.frame)
names(d.frame)[2] <- paste(sid)


			} else { d.frame <- "NULL" };

assign(paste("Nuclei", rj, sep=""), d.frame)

# for first value 
	if ( rj == 1 && d.frame != "NULL"  ) {
		first <- d.frame }
 
# combining successive data.frames
	if ( rj == 2 && d.frame != "NULL"  ) { 
		sum1 <- merge(first, d.frame, by="ID", all=TRUE) }

# combining successive data.frames
	if ( rj > 2 && d.frame != "NULL"  ) { 
		a <- count(1)
		sum1 <- merge(sum1, d.frame, by="ID", all=TRUE) 
		} else { print("skip") }

}


write.table(sum1,file=csvfile,sep=",",row.names=FALSE,col.names=TRUE)

#####################################
##### END OF READ VCF MODULE ########
#####################################

## Processing position information
## ie color label row side col as HCT, hTert or shared 

# combine into one file
hh<-rbind(hct116,htert)
all<-rbind(shared,hh)
######################################################################
# Read and plot the Temp Correlations for each file and out each as a PDF
for (i in 1:length(file_names)) {  

        # Filter out incomplete cases for flagging                          
        filt <- QC_list[[i]][!complete.cases(QC_list[[i]]),]
        # Include only complete cases with no NA values
        cases <- na.omit(QC_list[[i]])

        # Gets the sample type T,N,Xn etc from the name
        sample <- strsplit(colnames(QC_list[[i]]), split="_")[[2]][2]
        name <- paste(paste("Temp_Corr",SA,sep="_"), sample, sep="-")
        fname <- paste(path,name,sep="/")

        # Output each as a separate PDF
        pdf (file=paste(fname, ".pdf", sep=""))

        # Formats the plots
        title <- paste(paste("Correlation of Amplicon Melting Temperatures for", SA, sep=" "), sample, sep="-")
        plot(cases$Actual_Temp, cases$Cal_Temp, xlab="Actual Amplicon melting Temperature", ylab="Calculated Amplicon Melting Temp", main=title, pch=19, ylim=c(60,100), xlim=c(60,100))

        # This section prints out the wells that failed (missing data or no data)
        createCounter <- function(value) { function(j) { value <<- value+j} }
        counter <- createCounter(60)

                text(70,68, "Wells with missing data / no data:")
                for (j in 1:length(filt$PCR_well)) {x <- counter(round(40/length(filt$PCR_well))); text(x, 65, filt$PCR_well[j], cex=0.6); print(x)}

        # This labels every other point on the plot for QC purposes
        textxy(cases$Actual_Temp, cases$Cal_Temp, cases$PCR_well, cx=0.8, m=c(mean(cases$Actual_Temp),mean(cases$Cal_Temp)))

        dev.off() 

        print("end")

}

q()
