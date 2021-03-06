library(ROCR)
data <- read.table("ROC_input.txt",header=TRUE,row.names=1)
mymethod <- prediction(-data[["mymethod"]],data$truth)
ribotaper <- prediction(-data[["ribotaper"]],data$truth)
orfscore <- prediction(data[["ORFscore"]],data$truth)

mymethod1 = performance(mymethod,"tpr","fpr")
ribotaper1 = performance(ribotaper,"tpr","fpr")
orfscore1 = performance(orfscore,"tpr","fpr")

mymethod2 = performance(mymethod,"prec","rec")
ribotaper2 = performance(ribotaper,"prec","rec")
orfscore2 = performance(orfscore,"prec","rec")

pdf("ROC.pdf")
par(mfrow=c(2,1))
par(mai=c(1,2,0.1,2))
#
plot(mymethod1@x.values[[1]],mymethod1@y.values[[1]],xlab = mymethod1@x.name,ylab=mymethod1@y.name,col="red",type="l",ylim = c(0,1))
points(ribotaper1@x.values[[1]],ribotaper1@y.values[[1]],xlab = ribotaper1@x.name,ylab=ribotaper1@y.name,col="blue",type="l")
points(orfscore1@x.values[[1]],orfscore1@y.values[[1]],xlab = orfscore1@x.name,ylab=orfscore1@y.name,col="green",type="l")
legend(x="bottomright", legend=c("RiboCode", "RiboTaper", "ORFscore"),bty='n',col=c("red", "blue", "green"),lwd=2.5,cex=1.0)
#
plot(mymethod2@x.values[[1]],mymethod2@y.values[[1]],xlab = mymethod2@x.name,ylab=mymethod2@y.name,col="red",type="l",ylim = c(0,1)
)
points(ribotaper2@x.values[[1]],ribotaper2@y.values[[1]],xlab = ribotaper2@x.name,ylab=ribotaper2@y.name,col="blue",type="l")
points(orfscore2@x.values[[1]],orfscore2@y.values[[1]],xlab = orfscore2@x.name,ylab=orfscore2@y.name,col="green",type="l")
legend(x="bottomright", legend=c("RiboCode", "RiboTaper", "ORFscore"),bty='n',col=c("red", "blue", "green"),lwd=2.5,cex=1.0)
dev.off()
#
mymethod3 = performance(mymethod,"auc")
ribotaper3 = performance(ribotaper,"auc")
orfscore3 = performance(orfscore,"auc")

auc=c(mymethod3@y.values[[1]],ribotaper3@y.values[[1]],orfscore3@y.values[[1]])
write.table(auc,"auc.txt",sep="\t",row.names=F,col.names=F)
