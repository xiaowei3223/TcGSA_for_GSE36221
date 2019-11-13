################
# Get GEO data(GSE36221)
# AUTHOR: Chengshu Xie 
# Created: Sep.2019
# Last update: 
# GEO: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE36221
# PLATFORM: GPL6244
################

################################################
###### Convert probe IDs into gene symbols #####
################################################
library(hugene10sttranscriptcluster.db)
library(R.utils)
library(tidyverse)

get_input_data = function(series_matrix_file, pdata_file, geo_accession, bioc_annotation_package){  ### make sure of the R annotation packages
  
  if (file.exists(series_matrix_file) == TRUE & file.exists(pdata_file) == TRUE ) {
    cat("Reading available file...\n")
    exprSet = read.table(series_matrix_file, header=TRUE, row.names = 1, comment.char="!", sep='\t', quote="")
    pdata = read.table(pdata_file, header=TRUE, row.names = 1, comment.char="!", sep='\t', quote="")
  } else {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("GEOquery")
    library(GEOquery)
    geo_data = getGEO(geo_accession, destdir=".",getGPL = F)
    exprSet = exprs(geo_data[[1]])
    exprSet = as.data.frame(exprSet)
    pdata = pData(geo_data[[1]])    
  }
  
  # Convert probe IDs into gene symbols:
  exprSet$probe_id = rownames(exprSet)
  if (bioc_annotation_package == "hugene10sttranscriptcluster.db") { 
    probe2symbol_df = toTable(get("hugene10sttranscriptclusterSYMBOL"))}
  newmatrixdata = exprSet %>% 
    inner_join(probe2symbol_df,by="probe_id") %>% 
    select(-probe_id) %>%         
    select(symbol, everything()) %>% 
    mutate(rowMean = rowMeans(.[grep("GSM", names(.))])) %>% 
    arrange(desc(rowMean))  %>% 
    distinct(symbol,.keep_all = T) %>% 
    select(-rowMean) %>%
    tibble::column_to_rownames(colnames(.)[1])
  newmatrixdata = as.matrix(newmatrixdata)
  result = list("pdata" = pdata, "exprdata" = newmatrixdata)
  result
}

origin_data = get_input_data(series_matrix_file = "GSE36221_series_matrix.txt", 
                             pdata_file = "GSE36221_pdata.txt",
                             geo_accession = "GSE36221",
                             bioc_annotation_package = "hugene10sttranscriptcluster.db")


################################################
################# tidy data ####################
################################################
pdata = origin_data$pdata
pdata_de = subset(pdata, pdata[,40] != 3)
GSE36221 = origin_data$exprdata[, which(colnames(origin_data$exprdata) %in% rownames(pdata_de))]
GSE36221 = GSE36221[, order(colnames(GSE36221))]

pdata_mon0 = pdata_de[which(pdata_de[,39] == 0),]
pdata_mon6 = pdata_de[which(pdata_de[,39] == 6),]
pdata_mon30 = pdata_de[which(pdata_de[,39] == 30),]
data_new = list()
data_new$mon0 = t(subset(t(GSE36221), rownames(t(GSE36221)) %in% rownames(pdata_mon0)))
data_new$mon6 = t(subset(t(GSE36221), rownames(t(GSE36221)) %in% rownames(pdata_mon6)))
data_new$mon30 = t(subset(t(GSE36221), rownames(t(GSE36221)) %in% rownames(pdata_mon30)))

################################################
############### GSVA analysis ##################
################################################
library(GSVAdata)
library(GSEABase)
library(GSVA)
result.gsva = function(matrix){
  
  ###### Geneset file
  genesetcollection = getGmt("c2.cp.kegg.v7.0.symbols.gmt",
                             geneIdType = EntrezIdentifier(),
                             collectionType = BroadCollection(category="c2"))
  result.gsva = gsva(matrix, genesetcollection, method = "gsva", mx.diff = FALSE, parallel.sz=1, abs.ranking = FALSE, verbose=TRUE)
  result.plage = gsva(matrix, genesetcollection, method = "plage", mx.diff = FALSE, parallel.sz=1, abs.ranking = FALSE, verbose=TRUE)
  result.ssgsea = gsva(matrix, genesetcollection, method = "ssgsea", mx.diff = FALSE, parallel.sz=1, abs.ranking = FALSE, verbose=TRUE)
  result.zscore = gsva(matrix, genesetcollection, method = "zscore", mx.diff = FALSE, parallel.sz=1, abs.ranking = FALSE, verbose=TRUE)
  result = list("result.ssgsea" = result.ssgsea, "result.plage" = result.plage, 
                "result.gsva" = result.gsva, "result.zscore" = result.zscore)
}

result_gsva = lapply(data_new, result.gsva)

################################################
############ get p-value results ###############
################################################
library(limma)
run_limma = function(exprSet, pdata){
  f = factor(pdata[,40], levels=c("0", "1","2"))
  design = model.matrix(~0+f)
  colnames(design) = c("tre0", "tre1","tre2")
  rownames(design) = rownames(pdata) 
  contrast.matrix = makeContrasts(tre2-tre0, tre1-tre0, tre2-tre1,levels=design)
  fit = lmFit(exprSet, design)
  fit2 = contrasts.fit(fit, contrast.matrix)
  fit2 = eBayes(fit2)
  results = topTable(fit2, adjust="BH", n = Inf)
  results
}

p.value_gsva = list()
p.value_gsva_all = list()
for (i in 1:3) { 
  pdata_mon = switch(i,
                     pdata_mon0,
                     pdata_mon6,
                     pdata_mon30)
  for (j in 1:4) { 
    p.value_gsva[[j]] = run_limma(result_gsva[[i]][[j]], pdata_mon)
  }
  names(p.value_gsva) = c("ssgsea", "plage", "gsva", "zscore")
  p.value_gsva_all[[i]] = p.value_gsva
}
names(p.value_gsva_all) = c("mon0", "mon6", "mon30")


