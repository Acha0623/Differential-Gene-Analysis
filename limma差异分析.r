setwd('~/todo/data/早产/GSE73685')
#save(gse_1, file = "GSE73685.RData")
#load("GSE73685.RData")


library(GEOquery)
library(limma)
#library(tidyr)
#library(dplyr)
library(tidyverse)

# 1.1 GEO数据获取
## GSE24129 sample 8v8
options( 'download.file.method.GEOquery' = 'libcurl' )   #国区下载geo数据
gse_1 <- getGEO("GSE73685",                              #可能需要更改内存Sys.setenv("VROOM_CONNECTION_SIZE" = 131072*5)
                 file = "./GSE73685_series_matrix.txt.gz", #已经自行下载了   
                 #GSEMatrix = TRUE, #通过getGEO从网络下载
                 destdir = ".",getGPL = T,AnnotGPL = T)
                 
exprs_1 <- exprs(gse_1) #表达信息        #exprs(gse_1[[1]])有可能有多组的情况下
pdata_1 <- pData(gse_1) #样品处理信息
fdata_1 <- fData(gse_1) #平台注释信息

# 1.2 芯片数据前期处理
if(T){
  #最新版本差异表达数据前期处理
  fdata_raw <- data.frame(ID = fdata_1$ID, Gene_symbol = fdata_1$`Gene symbol`)
  #探针去除没有gene symbol的行
  fdata_raw <-fdata_raw[fdata_raw$Gene_symbol!="",]
  # 多个探针对应一个基因,所有///的元素都要
  con_probe_gene <- apply(fdata_raw, 1, function(x){          
    # 将探针号和基因号，进行---的链接
    paste(x[1],str_split(x[2],'///', simplify=T), sep = "---")
  })
  con_probe_gene <- data.frame(con = unlist(con_probe_gene))                 
  fdata_new <- separate(con_probe_gene,con, c("ID","Gene_symbol"),sep = '---')    
  dim(fdata_new) #28819
  #探针去重
  table(duplicated(fdata_new)) #633
  fdata_new <- fdata_new[!duplicated(fdata_new),] 
  dim(fdata_new)#28186
  #开始merge
  exprs_1 = data.frame(exprs_1)
  exprs_1$ID = rownames(exprs_1)
  mergeexprs_1 <- merge(exprs_1,fdata_new,by="ID")
  dim(mergeexprs_1) #28186
  # 去重(多个探针对应一个基因)
  table(duplicated(mergeexprs_1["Gene_symbol"])) #7039
  uniqexprs_1<- avereps(mergeexprs_1[,2:(dim(mergeexprs_1)[2]-1)],
                        ID=mergeexprs_1$`Gene_symbol`)
  dim(uniqexprs_1) #21147
}


## 分组
group_info_1 <- data.frame(cbind("sample"=pdata_1[,2],
                                 "tissue"= pdata_1$`tissue type:ch1`,
                                 "people" = pdata_1$`subject id:ch1`,
                                 "group_raw" = pdata_1$`outcome:ch1`))
group_info_1$group[grepl("preterm",group_info_1$group_raw)] = "Preterm"
group_info_1$group[grepl(" term",group_info_1$group_raw)] = "CTRL"
group_info <- group_info_1[group_info_1$group == "Preterm" |  group_info_1$group == "CTRL",]

## 组织分组 仅以Amnion为例
group_Amnion <- group_info[group_info$tissue == "Amnion",]


# 1.3 差异分析  Preterm-CTRL需要修改
deg_limma <- function(exprSet,group_list){
  # design
  design <- model.matrix(~0+factor(group_list))
  colnames(design)=levels(factor(group_list))
  rownames(design)=colnames(group_list)
  # contrast.matrix 
  contrast.matrix<-makeContrasts("Preterm-CTRL",levels = design)
  # fit model
  fit <- lmFit(exprSet,design)
  fit2 <- contrasts.fit(fit, contrast.matrix) 
  fit2 <- eBayes(fit2) 
  DEG<-topTable(fit2, coef=1, n=Inf) %>% na.omit() 
  DEG
}

## Amnion组织差异表达
uniqexprs_Amnion <- uniqexprs_1[,colnames(uniqexprs_1) %in% group_Amnion$sample] #20
deg_Amnion <- deg_limma(uniqexprs_Amnion,group_Amnion$group)
deg_sign_Amnion <- deg_Amnion[deg_Amnion$P.Value<0.05,] #961 

# 1.4 保存
write.csv(deg_sign_Amnion, "./DEG_sign_Amnion.csv", row.names = T)
write.table(deg_sign_Amnion, "./DEG_sign_Amnion.txt", sep = "\t",row.names = T, quote = F)


