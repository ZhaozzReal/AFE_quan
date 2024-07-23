suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(DEXSeq))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(plyr))
args = commandArgs(TRUE)
option.list = list(
  make_option(c("-b","--bam"),action = "store_true",default = FALSE,help = "input file contains all filenames of bamfile [%default]"),
  make_option(c("-c","--count"),action = "store_true",default = FALSE,help = "input file of FE expression count matrix generated in the previous step [%default]"),
  make_option(c("-d","--dir"),action = "store_true",default = FALSE,help = "input directory name contains all exoncount files for DEXSeq model [%default]"),
  make_option(c("-o","--output"),action = "store_true",default = FALSE,help = "output final results[%default]")
)
desc = "Infer significantly dysregulated alternative first exons between conditions using DEXSeq model"
parser = OptionParser(option_list = option.list, description = desc)
opt = parse_args(parser, args = args, positional_arguments = TRUE)

calc_cpm = function(expr_mat){
  expr_mat[is.na(expr_mat)] = 0
  norm_factor = colSums(expr_mat) + 1
  return (t(t(expr_mat)/norm_factor)*10^6) %>% as.data.frame()
}
Obtain_MajorPromoter_usage = function(a){
  a = a[order(as.numeric(a$median),decreasing = TRUE),]
  SYMBOL = a$SYMBOL[1]
  num_promoters = nrow(a)
  majorP = a$anno[1]
  majorP_median = a$median[1]
  majorPromoter_usage = apply(a %>% dplyr::select(contains(filenames)),2,function(x)round((x[1]/(sum(x)+0.1)),3))
  tmp = c(SYMBOL = SYMBOL,num_promoters=num_promoters,majorP=majorP,majorP_median=majorP_median,majorPromoter_usage)
  return(tmp)
}


cfg_file = read.table(opt$args[1])
filenames = apply(cfg_file,1,function(x)strsplit(x,"=") %>% unlist() %>% .[2] %>% strsplit(",") %>% unlist()) %>% as.character()
condition_length = apply(cfg_file,1,function(x)strsplit(x,"=") %>% unlist() %>% .[2] %>% strsplit(",") %>% unlist() %>% length())
condition_list = c(paste0("condition1_",c(1:condition_length[1])),paste0("condition2_",c(1:condition_length[2])))
count_table = fread(opt$args[2],header = T) %>% rowwise() %>% dplyr::mutate(first_exon_region=str_replace(first_exon_region,":","_")) %>% dplyr::mutate(anno=paste0(genename,":",first_exon_region))
for (i in 1:length(filenames)){
  filename = filenames[i]
  condition = condition_list[i]
  count_table %>% dplyr::select(anno,filenames[i]) %>% 
    write.table(paste0(opt$args[3],condition,"_count.txt"),col.names = F,sep = "\t",quote = FALSE,row.names = F)
}

sampleTable = data.frame(row.names = filenames,condition = c(rep("condition1",condition_length[1]),rep("condition2",condition_length[2])))
countFiles = list.files(opt$args[3], full.names=TRUE)
dxd = DEXSeqDataSetFromHTSeq(countFiles,sampleData = sampleTable,design = ~ sample + exon + condition:exon)
dxd = estimateSizeFactors( dxd )
dxd = estimateDispersions( dxd )
dxd = testForDEU( dxd )
dexresult = DEXSeqResults( dxd ) %>% as.data.frame() %>% na.omit() %>%  dplyr::select("groupID","featureID","padj") %>% rowwise() %>% dplyr::mutate(featureID=str_replace(featureID,"E",":")) %>% dplyr::mutate(majorP=paste0(groupID,featureID)) %>% dplyr::select(majorP,padj)


count_table_normalized = count_table %>% tibble::column_to_rownames(var="anno") %>% dplyr::select(filenames) %>% calc_cpm() %>% as.data.frame()
count_table_normalized$median = apply(count_table_normalized,1,median)
count_table_normalized = count_table_normalized %>% dplyr::filter(median>0.25) %>% tibble::rownames_to_column(var="anno")
count_table_normalized$SYMBOL = count_table_normalized$anno %>% lapply(function(x)strsplit(x,":",fixed=T) %>% unlist() %>% .[1]) %>% unlist()
count_table_normalized = ddply(count_table_normalized,.(SYMBOL),function(x){if(nrow(x)!=1)return(x)})
count_table_ratio = ddply(count_table_normalized,.(SYMBOL),Obtain_MajorPromoter_usage)
count_table_ratio$mean_condition1 = count_table_ratio %>% dplyr::select(filenames[1:condition_length[1]]) %>% apply(2,as.numeric) %>% rowMeans()
count_table_ratio$mean_condition2 = count_table_ratio %>% dplyr::select(filenames[(condition_length[1]+1):length(filenames)]) %>% apply(2,as.numeric) %>% rowMeans()
count_table_ratio$diff = count_table_ratio$mean_condition2 - count_table_ratio$mean_condition1
count_table_ratio = inner_join(count_table_ratio,dexresult)
final = count_table_ratio %>% dplyr::mutate(change=ifelse(padj<0.05&abs(diff)>0.05,ifelse(diff>0,"Up","Down"),"Nochange"))
write.table(final,opt$args[4],quote = FALSE,sep = "\t",row.names = FALSE)