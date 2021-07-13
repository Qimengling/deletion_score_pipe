library("Biostrings")
args=commandArgs(T)
info<-read.csv(args[1],header=T)
###repeat filter
###MR
MR_freq<-rep(0,dim(info)[1])
if (file.info("del_MR.tsv")$size >0 )
{
data<-read.table("del_MR.tsv",header = T)
index<-which(data$Spacer<=7 & data$Subset==1)
if (length(index)>0)
{
data<-data[index,]
for (k in 1:dim(data)[1])
{
temp<-strsplit(as.character(data$Sequence_name[k]),split = "")[[1]]
index1<-which(temp=="_")
index2<-which(temp=="(")
temp1<-as.numeric(paste(temp[(index1+1):(index2-1)],collapse = ""))
MR_freq[temp1]<-MR_freq[temp1]+1
}}}

###APR
APR_base<-read.table("APR_base.txt",header=F)
APR_base<-APR_base$x
APR_freq<-rep(0,dim(info)[1])
if (file.info("del_APR.tsv")$size >0)
{
data<-read.table("del_APR.tsv",header=T)
seq<-data$Sequence
for (k in 1:length(seq))
{
if (seq[k] %in% APR_base)
{
temp<-strsplit(as.character(data$Sequence_name[k]),split = "")[[1]]
index1<-which(temp=="_")
index2<-which(temp=="(")
temp1<-as.numeric(paste(temp[(index1+1):(index2-1)],collapse = ""))    
APR_freq[temp1]<-APR_freq[temp1]+1
}}}

repeat_type<-c("DR","IR","GQ","STR","Z")
repeat_freq<-MR_freq+APR_freq
for (i in 1:length(repeat_type))
{
temp<-sprintf("./del_%s.tsv",repeat_type[i])
if (file.info(temp)$size>0)
{    
data<-read.table(temp,header=T)
for (k in 1:dim(data)[1])
{
temp<-strsplit(as.character(data$Sequence_name[k]),split = "")[[1]]
index1<-which(temp=="_")
index2<-which(temp=="(")
temp1<-as.numeric(paste(temp[(index1+1):(index2-1)],collapse = ""))
repeat_freq[temp1]<-repeat_freq[temp1]+1
}}}
repeat_freq<-repeat_freq/2

a<-read.table("del_flank.fa")
a<-a$V1
a<-as.character(a[seq(2,length(a),by=2)])

gc_freq<-rep(0,dim(info)[1])
for (i in 1:dim(info)[1])
{
temp<-(strsplit(as.character(a[(i-1)*2+1]),split=""))[[1]]
gc_freq[i]<-gc_freq[i]+length(which(temp=="G" | temp=="C"))/length(temp)
temp<-(strsplit(as.character(a[i*2]),split=""))[[1]]
gc_freq[i]<-gc_freq[i]+length(which(temp=="G" | temp=="C"))/length(temp)
}
gc_freq<-gc_freq/2

motif<-read.table("./data/motif_info.txt")
motif<-as.character(motif$x)
motif_freq<-rep(0,dim(info)[1])
for (k in 1:dim(info)[1])
{
  for (m in 1:length(motif))
  {
    temp<-matchPattern(motif[m], DNAString(a[(k-1)*2+1]), fixed = c(subject = TRUE, pattern = FALSE))
    motif_freq[k]<-length(temp@ranges@start)+motif_freq[k]
    temp<-matchPattern(motif[m], DNAString(a[k*2]), fixed = c(subject = TRUE, pattern = FALSE))
    motif_freq[k]<-length(temp@ranges@start)+motif_freq[k]
  }
}
motif_freq<-motif_freq/2

repeat_fre_del<-read.csv("./data/repeat_freq_for_each_del_all.csv")
repeat_fre_del<-repeat_fre_del$x
motif_freq_for_each_del<-read.csv("./data/motif_freq_for_each_deletion.csv",row.names = 1)
motif_freq_for_each_del<-motif_freq_for_each_del$x
gc_for_each_del<-read.csv("./data/gc_freq_for_each_deletion.csv")
gc_for_each_del<-gc_for_each_del$x

score<-matrix(0,nrow=dim(info)[1],ncol=3)
for (i in 1:dim(info)[1])
{
score[i,1]<-(length(which(repeat_fre_del<=repeat_freq[i]))-
                 0.5*length(which(repeat_fre_del==repeat_freq[i])))/42098
score[i,2]<-(length(which(motif_freq_for_each_del<=motif_freq[i]))-
                 0.5*length(which(motif_freq_for_each_del==motif_freq[i])))/42098
score[i,3]<-(length(which(gc_for_each_del<=gc_freq[i]))-
                 0.5*length(which(gc_for_each_del==gc_freq[i])))/42098
}
result<-data.frame(info,score=rowSums(score))
write.csv(result,"score_result.csv",row.names=F)


