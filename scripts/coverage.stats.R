## #!/home/rebmmab/programs/R-3.2.5/bin Rscript

### #!/usr/bin/env Rscript

# Usage : coverage-plot.R <coveragefile> <genomelength>
args <- commandArgs(trailingOnly = TRUE) # capture command line arguments
# Put command line arguments into variables
coveragefile <- args[1]
localgenomelength <- args[2]

seqname <- gsub("^.*\\/","" ,coveragefile, perl=T)
#seqname <- gsub("\\..*$", "", seqname, perl=T)
seqname <- gsub("\\.depth$", "", seqname, perl=T)
seqpath <- gsub("\\.depth","",coveragefile)

mainpath <- gsub("variants_new\\/.+$","assembly_new\\/",seqpath)
totals <- localgenomelength

coverage = read.table(coveragefile,header=F,comment.char = "")  #Read in file

names(coverage) <- c("sample", "position","coverage")
meancov.over1 <- sum(coverage$coverage)/nrow(coverage)


# "depth" file does not report 0 coverage positions - reintroduce these here
allpositions <- data.frame(c(1: totals[[1]]),stringsAsFactors = F)
names(allpositions) <- "position"
coverage <- merge(coverage,allpositions, by="position", all.y=T)
coverage$coverage <- as.numeric(coverage$coverage)
coverage$coverage <- as.numeric(ifelse(is.na(coverage$coverage),"0",coverage$coverage))
coverage$sample <- seqname

genomelength = nrow(coverage)

# Generate and output coverage statistics
meancov.all <- round(sum(coverage$coverage)/nrow(coverage),1)
#meancov.all <- round(mean(coverage,na.rm=T),1)
mediancov.all <- median(coverage$coverage)
max.cov <- max(coverage$coverage)
min.cov <- min(coverage$coverage)
cov.250 <- round((sum(coverage$coverage >=250)/nrow(coverage))*100,1)
cov.100 <- round((sum(coverage$coverage >=100)/nrow(coverage))*100,1)
cov.10 <- round((sum(coverage$coverage >=10)/nrow(coverage))*100,1)
cov.8 <- round((sum(coverage$coverage >=8)/nrow(coverage))*100,1)
cov.5 <- round((sum(coverage$coverage >=5)/nrow(coverage))*100,1)
cov.4 <- round((sum(coverage$coverage >=4)/nrow(coverage))*100,1)
cov.2 <- round((sum(coverage$coverage >=2)/nrow(coverage))*100,1)
covstats <- cbind(seqname,genomelength,meancov.all,mediancov.all,max.cov,min.cov,cov.250,cov.100,cov.10,cov.8,cov.5,cov.4,cov.2)
write.table(covstats,file=paste(seqpath,".coverage.stats",sep=""),quote=F,sep="\t",col.names=T,row.names=F)

