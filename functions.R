annotate_gene <- function(data){
  #gene name annotation
  (gpl=data@annotation)
  checkGPL(gpl)
  printGPLInfo(gpl)
  probe2symbol_df <- idmap(gpl)
  length(unique(probe2symbol_df$probe_id))
  length(unique(probe2symbol_df$symbol))
  return(probe2symbol_df)
}

annotate_gene_from_db <-function(data){
  platformMap <- data.table::fread("./resources/platformMap.txt", data.table = F)
  index <- unique(pData(data)$platform)
  paste0(platformMap$bioc_package[grep(index,platformMap$gpl)],".db")
  paste("GEOID:",geoID, "Platform:",index)
  
  if(index == "GPL6244"){
    if(!requireNamespace("hugene10sttranscriptcluster.db")){
      options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
      BiocManager::install("hugene10sttranscriptcluster.db",update = F, ask = F)
    } 
    library(hugene10sttranscriptcluster.db)
    probe2symbol_df <- toTable(get("hugene10sttranscriptclusterSYMBOL"))
  }
  
  if(index == "GPL6947"){
    if(!requireNamespace("illuminaHumanv3.db")){
      options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
      BiocManager::install("hgu133a.db",update = F, ask = F)
    }
    library(illuminaHumanv3.db)
    probe2symbol_df <- toTable(get("illuminaHumanv3SYMBOL"))
  }
  
  if(index == "GPL96"){
    if(!requireNamespace("hgu133a.db")){
      options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
      BiocManager::install("hgu133a.db",update = F, ask = F)
    }
    library(hgu133a.db)
    probe2symbol_df <- toTable(get("hgu133aSYMBOL"))
  }
  if(index == "GPL5175"){
    if(!requireNamespace("huex10sttranscriptcluster.db")){
      options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
      BiocManager::install("huex10sttranscriptcluster.db",update = F, ask = F)
    }
    library(huex10sttranscriptcluster.db)
    probe2symbol_df <- toTable(get("huex10sttranscriptclusterSYMBOL"))
  }
  
  if(index == "GPL10558"){
    if(!requireNamespace("illuminaHumanv4.db")){
      options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
      BiocManager::install("illuminaHumanv4.db",update = F, ask = F)
    }
    library(illuminaHumanv4.db)
    probe2symbol_df <- toTable(get("illuminaHumanv4SYMBOL"))
  }
  
  if(index == "GPL571"){
    if(!requireNamespace("hgu133a2.db")){
      options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
      BiocManager::install("hgu133a2.db",update = F, ask = F)
    }
    library(hgu133a2.db)
    probe2symbol_df <- toTable(get("hgu133a2SYMBOL"))
  }
  
  length(unique(probe2symbol_df$probe_id))
  length(unique(probe2symbol_df$symbol))
  return(probe2symbol_df)
}
