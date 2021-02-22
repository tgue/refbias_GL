
library(ggplot2)

tab=read.table('q.summary.4X.txt')
tab$method=paste(tab$V8,tab$V7)

pdf('est.q.plot.4X.pdf',height=15,width=15)
ggplot(tab,aes(x=V1,y=V5,fill=method))+geom_boxplot()+ facet_grid(rows=vars(V2),cols=vars(V3))
dev.off()

