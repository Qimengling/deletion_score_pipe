###get deletion flanking bed file
args=commandArgs(T)
info<-read.csv(args[1],header=T)
bed_file<-matrix(nrow=dim(info)[1]*2,ncol=6)
for (i in 1:dim(info)[1])
{
bed_file[(i-1)*2+1,]<-c(paste0("chr",info[i,1]),
              as.numeric(as.character(info[i,2]))-500-1,
              as.numeric(as.character(info[i,2]))+499,
              paste("b1",i,sep="_"),
              1,
              as.character(info[i,4]))
bed_file[i*2,]<-c(paste0("chr",info[i,1]),
              as.numeric(as.character(info[i,3]))-499-1,
              as.numeric(as.character(info[i,2]))+500,
              paste("b2",i,sep="_"),
              1,
              as.character(info[i,4]))
}
write.table(bed_file,"bed_file.txt",row.names=F,col.names=F)
