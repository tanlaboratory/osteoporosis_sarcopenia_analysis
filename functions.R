annotate_gene <- function(data){
  #基因名注释
  (gpl=data@annotation)
  checkGPL(gpl)
  printGPLInfo(gpl)
  probe2symbol_df <- idmap(gpl)
  #probe2symbol_df
  
  ## 探针有多少个？
  length(unique(probe2symbol_df$probe_id))
  ## 这么多行中，基因名称有重复的么？
  length(unique(probe2symbol_df$symbol))
  return(probe2symbol_df)
}


annotate_gene_from_db <-function(data){
  #################################################################
  ## 探针基因名转换
  ##platformMap 中有常见的平台个R注释包的对应关系，这是我整理的。
  ## 读取，这都是我们已经讲过的
  platformMap <- data.table::fread("./resources/platformMap.txt", data.table = F)
  
  ## 平台的名称如何知道?
  index <- unique(pData(data)$platform)
  #index
  ## 数据储存在bioc_package这一列中
  paste0(platformMap$bioc_package[grep(index,platformMap$gpl)],".db")
  paste("GEOID:",geoID, "Platform:",index)
  
  ## 安装R包,可以直接安装，这里用了判断
  if(index == "GPL6244"){
    if(!requireNamespace("hugene10sttranscriptcluster.db")){
      options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
      BiocManager::install("hugene10sttranscriptcluster.db",update = F, ask = F)
    } 
    ## 加载R包
    library(hugene10sttranscriptcluster.db)
    ## 获取探针和基因的对应关系：这是探针注释的关键步骤
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
    #ref: https://plos.figshare.com/articles/dataset/Bioconductor_R_AnnotationData_packages_and_available_symbols_for_the_selected_series_integration_/6261251
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
  
  ## 探针有多少个？
  length(unique(probe2symbol_df$probe_id))
  ## 这么多行中，基因名称有重复的么？
  length(unique(probe2symbol_df$symbol))
  return(probe2symbol_df)
}

# annotate_gene_manually <- function(){
#   library(GEOquery)
#   library(Biobase)
#   library(data.table)
#   
#   # 获取平台ID
#   platform_id <- annotation(gset[[1]])
#   platform_id
#   # 从GEO网站下载平台注释文件 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi
#   # 请根据具体平台ID下载相应的注释文件
#   
#   if(platform_id == "GPL4133"){
#     annot_file <- "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL4nnn/GPL4133/suppl/GPL4133%5Fold%5Fannotations.txt.gz"
#   }elseif(platform_id == "GPL2700"){
#     annot_file <- "ftp://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL2nnn/GPL2700/annot/GPL2700.annot.gz"
#   }else{
#     annot_file <- paste0("ftp://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/", platform_id, "/annot/", platform_id, ".annot.gz")
#   }
#   annot_file
#   
#   dest_file  <- paste0("./data/", platform_id, "_annotation.gz")
#   dest_file
#   download.file(annot_file, destfile = dest_file)
#   
#   all_data <- readLines(dest_file)
#   head(all_data,30)
#   View(all_data)
#   # 使用grep找到数据起始行
#   data_start <- grep("ID\t", all_data)  # 假设数据部分以 "ID" 列开头
#   data_start <- data_start[1]
#   data_start
#   # 读取实际数据部分
#   #annot_data <- read.delim("annotation.gz", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)
#   annot_data <- fread(dest_file, skip = data_start-1, header = TRUE, sep = "\t", quote = "")
#   View(annot_data)
#   
#   # 假设注释文件中有一列“ID”表示探针ID，列“Gene.symbol”表示基因符号
#   # 合并表达数据和注释信息
#   id_column_index   <- which(colnames(annot_data) == "ID")
#   gene_symbol_index <- which(colnames(annot_data) %in% c("Gene symbol","Gene Symbol"))
#   gene_symbol_index
#   
#   probe2symbol_df <- annot_data[ ,c(id_column_index, gene_symbol_index)] # why this does not work?
#   probe2symbol_df <- annot_data[,c(1,3)]
#   colnames(probe2symbol_df) <- c("probe_id", "symbol")
#   
#   return(probe2symbol_df)
# }
