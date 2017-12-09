library(ROCR)
Args=commandArgs()
outname = Args[7]
data <- read.table(Args[6],header=FALSE,na.strings = "None",) 
if (outname == "ribocode"){
	# if ribocode: 
	pre  <- prediction(-data$V2,data$V3)
}else{
	pre  <- prediction(data$V2,data$V3) 
}
per <- performance(pre,"tpr","fpr")
prec <- performance(pre,"prec","rec")
pdf(paste0(outname,"_roc.pdf"))
plot(per)
dev.off()
auc <- performance(pre,"auc")
print (c(outname,auc@y.values[[1]]))

