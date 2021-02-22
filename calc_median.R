tab=read.table('q.summary.1X.txt')
tab$method=paste(tab$V8,tab$V7)

refs=unique(tab$V1)
prop=unique(tab$V2)
div=unique(tab$V3)
methods=unique(tab$method)



for (ref in refs){
s1=subset(tab,tab$V1==ref)
for (p in prop){
s2=subset(s1,s1$V2==p)
for (d in div){
s3=subset(s2,s2$V3==d)
for (m in methods){
s4=subset(s3,s3$method==m)
med=median(s4$V5)
cat(ref,p,d,m,med,'\n')
}
}
}
}

