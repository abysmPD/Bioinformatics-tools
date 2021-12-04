library(WGCNA)
library(flashClust)
options(stringsAsFactors=FALSE)
allowWGCNAThreads()

#######   #################    ################   #######    
#                     导入数据
#######   #################    ################   #######   

#-----Load expression data

exprMat <- "diff_FPKM.txt"
dataExpr <- read.table(exprMat, sep='\t', row.names=1, header=T, 
                     quote="", comment="", check.names=F)

dim(dataExpr)

#####数据筛选#######
## 筛选中位绝对偏差前25%的基因，至少MAD大于0.01
## 筛选后会降低运算量，也会失去部分信息
## 也可不做筛选，使MAD大于0即可
#m.mad <- apply(dataExpr,1,mad)
#dataExpr <- dataExpr[which(m.mad > 
#                 max(quantile(m.mad, probs=seq(0, 1, 0.85))[2],0.01)),]
                 
## 转换为样品在行，基因在列的矩阵
datExpr0 <- as.data.frame(t(dataExpr))

## 检测缺失值
 gsg = goodSamplesGenes(datExpr0, verbose = 3)

if (!gsg$allOK){
  #可选地，输出被删除的基因和样本名称:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", 
                     paste(names(datExpr0)[!gsg$goodGenes], collapse = ",")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", 
                     paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ",")));
  # 从数据中去除有问题基因和样本:
  datExpr0 = dataExpr[gsg$goodGenes, gsg$goodSamples]
}

nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)

dim(datExpr0)


#-----Load trait data
datTraits <- read.table("sample_info.txt", sep='\t', header=T, row.names=1,
                          check.names=FALSE, comment='',quote="")
                          
#######   #################    ################   #######    
#                 Call sample outliers
#######   #################    ################   #######   

#-----样品树状图和性状
A=adjacency(t(datExpr0),type="signed hybrid")

#-----计算整个网络连通性
k=as.numeric(apply(A,2,sum))-1

#-----标准化的连通性
Z.k=scale(k)
thresholdZ.k=-4 
outlierColor=ifelse(Z.k<thresholdZ.k,"red","black")
sampleTree = flashClust(as.dist(1-A), method = "complete")

#-----将性征转化成颜色
traitColors=data.frame(numbers2colors(datTraits,signed=T))
dimnames(traitColors)[[2]]=paste(names(datTraits))
datColors=data.frame(outlier=outlierColor,traitColors)

#-----画样本系统树
pdf('1_Sample_dendrogram_and_trait_heatmap.pdf')
plotDendroAndColors(sampleTree,groupLabels=names(datColors),
                    colors=datColors,main="Sample dendrogram and trait heatmap")
dev.off()

#-----删除异常样本(有必要删除时可以调用)
remove.samples= Z.k<thresholdZ.k | is.na(Z.k)
datExpr0=datExpr0[!remove.samples,]
datTraits=datTraits[!remove.samples,]
A=adjacency(t(datExpr0),type="distance")
k=as.numeric(apply(A,2,sum))-1
Z.k=scale(k)

save(datExpr0, datTraits, file="SamplesAndTraits.RData")

#######   #################    ################   #######    
#                     选择软阈值
#######   #################    ################   #######     
options(stringsAsFactors = FALSE)
lnames= load(file="SamplesAndTraits.RData")
lnames
dim(datExpr0)
dim(datTraits)
powers= c(seq(1,10,by=1), seq(from =12, to=40, by=2)) #choosing a set of soft-thresholding powers
sft = pickSoftThreshold(datExpr0, powerVector=powers, verbose =5,networkType="signed hybrid") #call network topology analysis function

pdf('2_Soft_Threshold_select.pdf', width = 12, height=8)
par(mfrow= c(1,2))
cex1=0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab= "Soft Threshold (power)", 
                       ylab="Scale Free Topology Model Fit, signed hybrid", type= "n", main= paste("Scale independence"))

text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers, cex=cex1, col="red")
abline(h=0.9, col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab= "Soft Threshold (power)", ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1, col="red")
dev.off()

softPower=sft$powerEstimate           #还是要自己看一下   软阈值太低会使模块数量下降，单模块基因数上升
softPower = 14                
#######   #################    ################   #######    
#                    构建网络
#######   #################    ################   #######     

adjacency=adjacency(datExpr0, power=softPower, type="signed hybrid") 
TOM= TOMsimilarity(adjacency, TOMType="signed")     # TOMType:unsigned, signed, signed Nowick, unsigned 2, signed 2, signed Nowick 2
dissTOM= 1-TOM

geneTree= flashClust(as.dist(dissTOM), method="average")

sizeGrWindow(9, 5)
plot(geneTree, xlab="", sub="", main= "Gene Clustering on TOM-based dissimilarity", labels= FALSE, hang=0.04)

#######   #################    ################   #######    
#                    Make modules
#######   #################    ################   ####### 

minModuleSize=30 
dynamicMods= cutreeDynamic(dendro= geneTree, distM= dissTOM, deepSplit=2, pamRespectsDendro= FALSE, minClusterSize= minModuleSize,method = "hybrid")
table(dynamicMods)

dynamicColors= labels2colors(dynamicMods)


plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang= 0.05, main= "Gene dendrogram and module colors")

#-----Merge modules whose expression profiles are very similar
MEList= moduleEigengenes(datExpr0, colors= dynamicColors)
MEs= MEList$eigengenes
#Calculate dissimilarity of module eigenegenes
MEDiss= 1-cor(MEs)
#Cluster module eigengenes
METree= flashClust(as.dist(MEDiss), method= "average")


plot(METree, main= "Clustering of module eigengenes", xlab= "", sub= "")
MEDissThres= 0.2
abline(h=MEDissThres, col="red")
merge= mergeCloseModules(datExpr0, dynamicColors, cutHeight= MEDissThres, verbose =3)

mergedColors= merge$colors
mergedMEs= merge$newMEs

pdf('3_merge_modules_tree.pdf', width = 12, height = 8)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang=0.05)
dev.off()

moduleColors= mergedColors
colorOrder= c("grey", standardColors(50))
moduleLabels= match(moduleColors, colorOrder)-1
MEs=mergedMEs

save(MEs, moduleLabels, moduleColors, geneTree, file= "SamplesAndColors_thresh24merge42_signed.RData")

#模块间的相关性
sizeGrWindow(12, 12) 
plotEigengeneNetworks(MEs, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2), plotDendrograms = FALSE, xLabelsAngle = 90)
#plotMEpairs(MEs,y=datTraits$cellType) #一般不运行
#######   #################    ################   #######    
#                Relate modules to traits
#######   #################    ################   ####### 

datt=datExpr0

#-----Define numbers of genes and samples
nGenes = ncol(datt);
nSamples = nrow(datt);
#-----Recalculate MEs with color labels
MEs0 = moduleEigengenes(datt, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

#-----Correlations of genes with eigengenes
moduleGeneCor=cor(MEs,datt)
moduleGenePvalue = corPvalueStudent(moduleGeneCor, nSamples);

moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
#---------------------Module-trait heatmap


textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
pdf('4_module_correlation_heatmap.pdf', width = 10, height = 15)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.7,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()
######--------------------end--------------------#######



#找到模块后如何筛选hub gene？
#1 High intramodular k within the module(KIM)
#2 High module membership (kMM,表达值与ME高相关)

# 计算模块成员值 (aka. 基于模块特征基因的连接kME):
datKME = signedKME(datExpr0, MEs)
#等同于
MM= as.data.frame(cor(datExpr0, MEs, use ="p"))

#计算模块内部gene连接度
KIM = intramodularConnectivity(adjacency, moduleColors, scaleByMax= TRUE)



#---------------------Gene significance by Module membership scatterplots
whichTrait="TP_type" #Replace this with the trait of interest


nGenes = ncol(datt);
nSamples = nrow(datt);
selTrait = as.data.frame(datTraits[,whichTrait]);
names(selTrait) = whichTrait
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(signedKME(datt, MEs));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datt, selTrait, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(selTrait), sep="");
names(GSPvalue) = paste("p.GS.", names(selTrait), sep="");
par(mfrow=c(2,3))
counter=0
for(module in modNames[1:length(modNames)]){
  counter=counter+1
  if (counter>6) {
    quartz()
    par(mfrow=c(2,3))
    counter=1
  }
  column = match(module, modNames);
  moduleGenes = moduleColors==module;
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste(module,"module membership"),
                     ylab = paste("GS for", whichTrait),
                     col = module,mgp=c(2.3,1,0))
}
######-------------------end------------------#######
write.table(merge[["colors"]],"color2.txt",sep = "\t")
write.table(datKME,"DATTKIM2.txt",sep = "\t")



#画单模块表达热图
#---------------------Eigengene heatmap
which.module="coral1" #replace with module of interest
datME=MEs
datExpr=datt

dir.create("MEs_heatmap")

ME=datME[, paste("ME",which.module, sep="")]
pdf(paste("MEs_heatmap/",which.module, '_ME.pdf', sep = ''), height = 5, width = 6.5)
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(datExpr[,moduleColors==which.module ]) ),
        nrgcols=30,rlabels=F,rcols=which.module,
        main=which.module, cex.main=2)
par(mar=c(5, 4.2, 0, 0.7))
barplot(ME, col=which.module, main="", names.arg=c(row.names(datt)), cex.names=0.5, cex.main=2,
        ylab="eigengene expression",xlab="sample")
dev.off()
######--------------------end--------------------#######
#画所有模块表达热图
#---------------------Eigengene heatmap
all.module=unique(mergedColors)

dir.create("MEs_heatmap")

for (m in all.module) {
  which.module=m #replace with module of interest
  datME=MEs
  datExpr=datt
  
  ME=datME[, paste("ME",which.module, sep="")]
  pdf(paste("MEs_heatmap/",which.module, '_ME.pdf', sep = ''), height = 5, width = 7)
  par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
  plotMat(t(scale(datExpr[,moduleColors==which.module ]) ),
          nrgcols=30,rlabels=F,rcols=which.module,
          main=which.module, cex.main=2)
  par(mar=c(5, 4.2, 0, 0.7))
  barplot(ME, col=which.module, main="", names.arg=c(row.names(datt)), cex.names=0.5, cex.main=2,
          ylab="eigengene expression",xlab="sample")
  dev.off()
}

######--------------------end--------------------#######



#######   #################    ################   #######    
#             导出到cytoscape 的信息
#######   #################    ################   ####### 
modules = c("grey60")
inModule=is.finite(match(moduleColors,modules))
probes = rownames(dataExpr)        #提取基因名
modProbes=probes[inModule]

# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files for Cytoscape
cyt = exportNetworkToCytoscape(modTOM,
  edgeFile=paste("CytoEdge",paste(modules,collapse="-"),".txt",sep=""),
  nodeFile=paste("CytoNode",paste(modules,collapse="-"),".txt",sep=""),
  weighted = TRUE, threshold = 0,nodeNames=modProbes,
  #altNodeNames = modGenes, 
  nodeAttr = moduleColors[inModule])

#######   #################    ################   #######    
#             批量导出到cytoscape 的信息
#######   #################    ################   ####### 
all.module=unique(mergedColors)

dir.create("MEs_edge_node")

for (m in all.module) {
modules = c(m)
inModule=is.finite(match(moduleColors,modules))
probes = rownames(dataExpr)        #提取基因名
modProbes=probes[inModule]

# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files for Cytoscape
cyt = exportNetworkToCytoscape(modTOM,
  edgeFile=paste("MEs_edge_node/CytoEdge",paste(modules,collapse="-"),".txt",sep=""),
  nodeFile=paste("MEs_edge_node/CytoNode",paste(modules,collapse="-"),".txt",sep=""),
  weighted = TRUE, threshold = 0,nodeNames=modProbes,
  #altNodeNames = modGenes, 
  nodeAttr = moduleColors[inModule])
}
































#######   #################    ################   #######    
#             Gene expression within modules
#######   #################    ################   ####### 



#---------------------Heatmap for top-kME in a module with gene names
vsd=read.csv("~/Desktop/diseaseScript_Final/VSDandPVALS_disease.csv")
names(vsd)

a.vsd=vsd[c(2:25)] #Columns with vsd
row.names(a.vsd)=vsd$X
head(a.vsd)
a.vsd=a.vsd[c(1:20,22:24)] #Remove wgcna outlier
names(a.vsd)

allkME =as.data.frame(signedKME(t(a.vsd), MEs))

gg=read.table("~/Documents/genomes/ahya_annotations_may25_2014/ahya2digNvec_plus_iso2gene.tab", sep="\t")
head(gg)
library(pheatmap)

whichModule="green"
top=25

modcol=paste("kME",whichModule,sep="")
sorted=a.vsd[order(allkME[,modcol],decreasing=T),]
hubs=sorted[1:top,]
# attaching gene names
gnames=c();counts=0
for(i in 1:length(hubs[,1])) {
  if (row.names(hubs)[i] %in% gg$V1) { 
    counts=counts+1
    gn=gg[gg$V1==row.names(hubs)[i],2]
    if (gn %in% gnames) {
      gn=paste(gn,counts,sep=".")
    }
    gnames=append(gnames,gn) 
  } else { 
    gnames=append(gnames,i)
  }
} 
row.names(hubs)=gnames

contrasting = colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)
contrasting2 = colorRampPalette(rev(c("chocolate1","chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)
contrasting3 = colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan","cyan")))(100)


pheatmap(hubs,scale="row",col=contrasting,border_color=NA, main=paste(whichModule,"top",top,"kME",sep=""))
######--------------------end--------------------#######


