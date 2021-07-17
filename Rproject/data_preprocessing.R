mkdir<-function(d){
  if (! file.exists(d)){
    dir.create(d)
  }
}
cancer='bcc'
mkdir(paste0("../data/",cancer,"_scRNA-seq"))
setwd(paste0("../data/",cancer,"_scRNA-seq"))
### data download

#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123813

### data loading
# Count data
library(dplyr)
bcc_matrix=gzfile('GSE123813_bcc_scRNA_counts.txt.gz','rt')  
dat=read.csv(bcc_matrix,header=T,sep='\t')

# meta data
GSE123813_bcc_all_metadata <- read.delim("GSE123813_bcc_all_metadata.txt")
GSE123813_bcc_scRNA_cells <- colnames(dat)
sample_bcc=GSE123813_bcc_all_metadata[GSE123813_bcc_all_metadata[,'cell.id'] %in% GSE123813_bcc_scRNA_cells,]

### identify T cells by annotation
CD4<-c("CD4_T_cells")
sample_cd4t=sample_bcc[sample_bcc[,'cluster'] %in% CD4,]

CD8<-c("CD8_mem_T_cells","CD8_ex_T_cells","CD8_act_T_cells")
sample_cd8=sample_bcc[sample_bcc[,'cluster'] %in% CD8,]

### identify CD4 T cells by CD3/CD4 expression
dat_cd4=dat[,dat['CD3D',]>0 | dat['CD3E',]>0]
dat_cd4=dat_cd4[,dat_cd4['CD4',]>0]                   

### identify CD8 T cells by CD3/CD8 expression levels
dat_cd8=dat[,dat['CD8A',]>0 | dat['CD8B',]>0]
dat_cd8=dat_cd8[,dat_cd8['CD3D',]>0|dat_cd8['CD3E',]>0]                

### data normalization at the cell level
#   by total count
#   this step is to remove the batch effects of cells
#   We used the same strategy as Seurat (https://satijalab.org/seurat/) 
#   and Scanpy (https://scanpy.readthedocs.io/en/stable/)

data.countExpressed <-
  apply(dat_cd8, 2, function(x){sum(as.numeric(x[x>0]))})

dat_cd8 <- (t(t(dat_cd8) / data.countExpressed))*10000

### Data for the functional signature and variance for the markers
markers<-c('exhaustion','eff','mem')
for (marker in markers){
  signature <- read.csv(paste0("../signatures/",marker,".csv"), header=TRUE, skip=1)
  
  genes<-signature[,'Gene.Symbol']
  dat_cd8_marker<-dat_cd8[genes,]
  
  SD=apply(dat_cd8_marker,1, sd, na.rm = TRUE)
  marker_sd<-data.frame(genes,SD)
  out<-"../Marker_dat/"
  mkdir(out)
  out<-paste0(out,'/',cancer,'_',marker)
  mkdir(out)
  write.table(marker_sd,paste0(out,'/','allmarkers_var.txt'),row.names = FALSE,quote=FALSE,sep='\t')
  write.table(t(dat_cd8_marker),paste0(out,'/','heatmap_data.txt.',cancer),row.names = TRUE,quote=FALSE,sep='\t')
}