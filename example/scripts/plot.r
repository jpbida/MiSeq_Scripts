####################################################################################################
### plot.r 
###
### 2012 Jp Bida
###
### This script is really specific to the Add4 results and will need to be customized 
### for different types of experiments 
###
### inputs: output from makeAln.cpp
### outputs: images in ./results 
####################################################################################################


data<-read.table(file="../data/add4lib.dat",sep=",")
### Alignment Quality Stats ####
names(data)=c("eid","sid","pos","qual")
data$eid=as.character(data$eid)
data$sid=as.character(data$sid)
data$qual=as.numeric(data$qual)
eids=table(data$eid)
eids=sort(eids,decreasing=TRUE)

### Show Top 10 Categories and Group everything else ###
n=as.character(names(eids[1:9]))
v=as.numeric(eids[1:9])
vo=as.numeric(sum(eids[10:length(eids)]))
no="others"
eid_table=cbind(c(n,no),c(v,vo))

### Show distribution of alignment quality to 3' end ###
jpeg(file="../results/AlignmentQuality.jpeg",quality=100)
hist(data$qual,nclass=10,col=2,main="Alignment Quality",xlab="Edit Distance(0=perfect,-10=bad)")
dev.off()

### Set a threshold for quality and keep only mapped eid's###
k=data$qual > -4 & (data$eid=="NOMOD" | data$eid=="A1M7" | data$eid=="NA1M7_1" | data$eid=="NA1M7_2" | data$eid=="NA1M7_3")
clean=data[k,]


### Get Single Point Mutation Data ###
len=nchar(clean$sid)
singles=len < 6 & len > 3
ptmuts=clean[singles,]
### Get positions ###
pos=gsub("[AGTC:]","",ptmuts$sid,perl=TRUE)
pos=as.numeric(pos)
ptmuts$mpos=pos
jpeg(file="../results/MutantDistribution.jpeg",quality=100)
### Look at counts for each position ###
par(mfrow=c(5,1),mai=c(0,0,0.2,0))
hist(ptmuts$mpos[ptmuts$eid=="NOMOD"],nclass=max(pos),col=2,main="NOMOD")
hist(ptmuts$mpos[ptmuts$eid=="A1M7"],nclass=max(pos),col=2,main="A1M7")
hist(ptmuts$mpos[ptmuts$eid=="NA1M7_1"],nclass=max(pos),col=2,main="NA1M7_1")
hist(ptmuts$mpos[ptmuts$eid=="NA1M7_2"],nclass=max(pos),col=2,main="NA1M7_2")
hist(ptmuts$mpos[ptmuts$eid=="NA1M7_3"],nclass=max(pos),col=2,main="NA1M7_3")
dev.off()

jpeg(file="../results/A1M7.jpeg",quality=100)
bgscale=5
adat=ptmuts[ptmuts$eid=="A1M7",]
adat$pos[adat$pos < 0]=0
maxpos=max(ptmuts$mpos)
minpos=0
plot(c(0,(10*maxpos)),c(0,(10*maxpos)),type="n")
for(i in min(ptmuts$mpos):max(ptmuts$mpos)){
tdat=adat[adat$mpos==i,]
### Normalize to the total number of hits ##
counts=as.numeric(table(tdat$pos))
pos=as.numeric(names(table(tdat$pos)))
full_ext=length(tdat$pos[tdat$pos==0])
ncounts=NULL
ncols=NULL
for(j in c(0:maxpos)){
ncounts[j+1]=0
if(length(counts[pos==j]) > 0){
ncounts[j+1]=counts[pos==j]/(full_ext+sum(counts[pos>=j]))
}
}
scale=100/max(ncounts)
ncols=colors()[253-floor(pmin(100,ncounts*scale*bgscale))]
rect(seq(0,(length(ncols)*10),length=length(ncols)),rep((10*i),length(ncols)),seq(10,((length(ncols)+1)*10),length=length(ncols)),rep((10*(i+1)),length(ncols)),col=ncols,border=ncols)
}
dev.off()

jpeg(file="../results/NOMOD.jpeg",quality=100)
bgscale=10
adat=ptmuts[ptmuts$eid=="NOMOD",]
adat$pos[adat$pos < 0]=0
maxpos=max(ptmuts$mpos)
minpos=0
plot(c(0,(10*maxpos)),c(0,(10*maxpos)),type="n")
for(i in min(ptmuts$mpos):max(ptmuts$mpos)){
tdat=adat[adat$mpos==i,]
### Normalize to the total number of hits ##
counts=as.numeric(table(tdat$pos))
pos=as.numeric(names(table(tdat$pos)))
full_ext=length(tdat$pos[tdat$pos==0])
ncounts=NULL
ncols=NULL
for(j in c(0:maxpos)){
ncounts[j+1]=0
if(length(counts[pos==j]) > 0){
ncounts[j+1]=counts[pos==j]/(full_ext+sum(counts[pos>=j]))
}
}
scale=100/max(ncounts)
ncols=colors()[253-floor(pmin(100,ncounts*scale*bgscale))]
rect(seq(0,(length(ncols)*10),length=length(ncols)),rep((10*i),length(ncols)),seq(10,((length(ncols)+1)*10),length=length(ncols)),rep((10*(i+1)),length(ncols)),col=ncols,border=ncols)
}
dev.off()

jpeg(file="../results/NA1M7_3.jpeg",quality=100)
bgscale=3
adat=ptmuts[ptmuts$eid=="NA1M7_3",]
adat$pos[adat$pos < 0]=0
maxpos=max(ptmuts$mpos)
minpos=0
plot(c(0,(10*maxpos)),c(0,(10*maxpos)),type="n")
for(i in min(ptmuts$mpos):max(ptmuts$mpos)){
tdat=adat[adat$mpos==i,]
### Normalize to the total number of hits ##
counts=as.numeric(table(tdat$pos))
pos=as.numeric(names(table(tdat$pos)))
full_ext=length(tdat$pos[tdat$pos==0])
ncounts=NULL
ncols=NULL
for(j in c(0:maxpos)){
ncounts[j+1]=0
if(length(counts[pos==j]) > 0){
ncounts[j+1]=counts[pos==j]/(full_ext+sum(counts[pos>=j]))
}
}
scale=100/max(ncounts)
ncols=colors()[253-floor(pmin(100,ncounts*scale*bgscale))]
rect(seq(0,(length(ncols)*10),length=length(ncols)),rep((10*i),length(ncols)),seq(10,((length(ncols)+1)*10),length=length(ncols)),rep((10*(i+1)),length(ncols)),col=ncols,border=ncols)
}
dev.off()

jpeg(file="../results/NA1M7_2.jpeg",quality=100)
bgscale=3
adat=ptmuts[ptmuts$eid=="NA1M7_2",]
adat$pos[adat$pos < 0]=0
maxpos=max(ptmuts$mpos)
minpos=0
plot(c(0,(10*maxpos)),c(0,(10*maxpos)),type="n")
for(i in min(ptmuts$mpos):max(ptmuts$mpos)){
tdat=adat[adat$mpos==i,]
### Normalize to the total number of hits ##
counts=as.numeric(table(tdat$pos))
pos=as.numeric(names(table(tdat$pos)))
full_ext=length(tdat$pos[tdat$pos==0])
ncounts=NULL
ncols=NULL
for(j in c(0:maxpos)){
ncounts[j+1]=0
if(length(counts[pos==j]) > 0){
ncounts[j+1]=counts[pos==j]/(full_ext+sum(counts[pos>=j]))
}
}
scale=100/max(ncounts)
ncols=colors()[253-floor(pmin(100,ncounts*scale*bgscale))]
rect(seq(0,(length(ncols)*10),length=length(ncols)),rep((10*i),length(ncols)),seq(10,((length(ncols)+1)*10),length=length(ncols)),rep((10*(i+1)),length(ncols)),col=ncols,border=ncols)
}
dev.off()

jpeg(file="../results/NA1M7_1.jpeg",quality=100)
bgscale=3
adat=ptmuts[ptmuts$eid=="NA1M7_1",]
adat$pos[adat$pos < 0]=0
maxpos=max(ptmuts$mpos)
minpos=0
plot(c(0,(10*maxpos)),c(0,(10*maxpos)),type="n")
for(i in min(ptmuts$mpos):max(ptmuts$mpos)){
tdat=adat[adat$mpos==i,]
### Normalize to the total number of hits ##
counts=as.numeric(table(tdat$pos))
pos=as.numeric(names(table(tdat$pos)))
full_ext=length(tdat$pos[tdat$pos==0])
ncounts=NULL
ncols=NULL
for(j in c(0:maxpos)){
ncounts[j+1]=0
if(length(counts[pos==j]) > 0){
ncounts[j+1]=counts[pos==j]/(full_ext+sum(counts[pos>=j]))
}
}
scale=100/max(ncounts)
ncols=colors()[253-floor(pmin(100,ncounts*scale*bgscale))]
rect(seq(0,(length(ncols)*10),length=length(ncols)),rep((10*i),length(ncols)),seq(10,((length(ncols)+1)*10),length=length(ncols)),rep((10*(i+1)),length(ncols)),col=ncols,border=ncols)
}
dev.off()
