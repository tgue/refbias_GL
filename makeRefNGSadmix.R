args<-commandArgs(trailingOnly = T)

if(length(args)==0){

    print("Arguments have to be supplied: ")
    print("1. first 3 cols of beagle file, 2. .filter file from NGSadmix")
    q()
}

tmpFile<-args[1]
filterFile<-args[2]

tmp<-read.table(tmpFile,as.is=T,h=T)
filter<-read.table(filterFile,as.is=T,h=T)

tmp2<-tmp[ tmp[,1]%in%filter[,1],]

ref<-cbind(id=tmp2[,1],chr=sapply(tmp2[,1],function(x) unlist(strsplit(x,"_"))[1]),pos=sapply(tmp2[,1],function(x) unlist(strsplit(x,"_"))[2]),name=tmp2[,1],A0_freq=tmp2[,3],A1=tmp2[,2])

write.table(ref,"tmp.ref",col=F,row=F,qu=F)
