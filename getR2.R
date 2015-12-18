args <-commandArgs()
myfile<-args[6]

correlaciones<-read.table(myfile,header=FALSE)

summary(lm(V1~V2, data=correlaciones))
