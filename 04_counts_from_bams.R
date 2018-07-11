# File: 04_counts_from_bams.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: generate count tables for transcripts from bam files
# Date: 11/07/2018


## set variables and source libraries
source('header.R')

## load the transcript db objects
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(GenomicAlignments)
library(rtracklayer)
# get the exons into GRangesList object
oGRLgenes = exonsBy(TxDb.Mmusculus.UCSC.mm10.knownGene, by = 'gene')

## the chromosome names used by the sequencing center in their index
# are not concordant with the standard naming, hence we rename all the genes
# in the oGRLgenes object

# > table(seqnames(unlist(oGRLgenes)))
# 
# chr1                 chr2                 chr3                 chr4                 chr5 
# 15404                19888                11267                14985                15690 
# chr6                 chr7                 chr8                 chr9                chr10 
# 12082                18184                11833                13924                11371 
# chr11                chr12                chr13                chr14                chr15 
# 18513                 8216                 8263                 8801                 9765 
# chr16                chr17                chr18                chr19                 chrX 
# 7367                11620                 5983                 7973                 8902 
# chrY                 chrM chr1_GL456210_random chr1_GL456211_random chr1_GL456212_random 
# 641                    0                    0                   13                    0 
# chr1_GL456213_random chr1_GL456221_random chr4_GL456216_random chr4_GL456350_random chr4_JH584292_random 
# 0                   10                    7                   21                    6 
# chr4_JH584293_random chr4_JH584294_random chr4_JH584295_random chr5_GL456354_random chr5_JH584296_random 
# 35                   28                    0                   17                    5 
# chr5_JH584297_random chr5_JH584298_random chr5_JH584299_random chr7_GL456219_random chrX_GL456233_random 
# 5                    8                   22                    4                   20 
# chrY_JH584300_random chrY_JH584301_random chrY_JH584302_random chrY_JH584303_random       chrUn_GL456239 
# 0                    0                    0                    0                    0 
# chrUn_GL456359       chrUn_GL456360       chrUn_GL456366       chrUn_GL456367       chrUn_GL456368 
# 0                    0                    0                    0                    0 
# chrUn_GL456370       chrUn_GL456372       chrUn_GL456378       chrUn_GL456379       chrUn_GL456381 
# 0                    0                    0                    0                    0 
# chrUn_GL456382       chrUn_GL456383       chrUn_GL456385       chrUn_GL456387       chrUn_GL456389 
# 0                    0                    0                    0                    0 
# chrUn_GL456390       chrUn_GL456392       chrUn_GL456393       chrUn_GL456394       chrUn_GL456396 
# 0                    0                    0                    0                    0 
# chrUn_JH584304 
# 11 
# > bf = BamFile(file.choose())
# > seqnames(seqinfo(bf))
# [1] "1"          "10"         "11"         "12"         "13"         "14"         "15"         "16"        
# [9] "17"         "18"         "19"         "2"          "3"          "4"          "5"          "6"         
# [17] "7"          "8"          "9"          "MT"         "X"          "Y"          "JH584299.1" "GL456233.1"
# [25] "JH584301.1" "GL456211.1" "GL456350.1" "JH584293.1" "GL456221.1" "JH584297.1" "JH584296.1" "GL456354.1"
# [33] "JH584294.1" "JH584298.1" "JH584300.1" "GL456219.1" "GL456210.1" "JH584303.1" "JH584302.1" "GL456212.1"
# [41] "JH584304.1" "GL456379.1" "GL456216.1" "GL456393.1" "GL456366.1" "GL456367.1" "GL456239.1" "GL456213.1"
# [49] "GL456383.1" "GL456385.1" "GL456360.1" "GL456378.1" "GL456389.1" "GL456372.1" "GL456370.1" "GL456381.1"
# [57] "GL456387.1" "GL456390.1" "GL456394.1" "GL456392.1" "GL456382.1" "GL456359.1" "GL456396.1" "GL456368.1"
# [65] "JH584292.1" "JH584295.1" "ERCC-00002" "ERCC-00003" "ERCC-00004" "ERCC-00009" "ERCC-00012" "ERCC-00013"
# [73] "ERCC-00014" "ERCC-00016" "ERCC-00017" "ERCC-00019" "ERCC-00022" "ERCC-00024" "ERCC-00025" "ERCC-00028"
# [81] "ERCC-00031" "ERCC-00033" "ERCC-00034" "ERCC-00035" "ERCC-00039" "ERCC-00040" "ERCC-00041" "ERCC-00042"
# [89] "ERCC-00043" "ERCC-00044" "ERCC-00046" "ERCC-00048" "ERCC-00051" "ERCC-00053" "ERCC-00054" "ERCC-00057"
# [97] "ERCC-00058" "ERCC-00059" "ERCC-00060" "ERCC-00061" "ERCC-00062" "ERCC-00067" "ERCC-00069" "ERCC-00071"
# [105] "ERCC-00073" "ERCC-00074" "ERCC-00075" "ERCC-00076" "ERCC-00077" "ERCC-00078" "ERCC-00079" "ERCC-00081"
# [113] "ERCC-00083" "ERCC-00084" "ERCC-00085" "ERCC-00086" "ERCC-00092" "ERCC-00095" "ERCC-00096" "ERCC-00097"
# [121] "ERCC-00098" "ERCC-00099" "ERCC-00104" "ERCC-00108" "ERCC-00109" "ERCC-00111" "ERCC-00112" "ERCC-00113"
# [129] "ERCC-00116" "ERCC-00117" "ERCC-00120" "ERCC-00123" "ERCC-00126" "ERCC-00130" "ERCC-00131" "ERCC-00134"
# [137] "ERCC-00136" "ERCC-00137" "ERCC-00138" "ERCC-00142" "ERCC-00143" "ERCC-00144" "ERCC-00145" "ERCC-00147"
# [145] "ERCC-00148" "ERCC-00150" "ERCC-00154" "ERCC-00156" "ERCC-00157" "ERCC-00158" "ERCC-00160" "ERCC-00162"
# [153] "ERCC-00163" "ERCC-00164" "ERCC-00165" "ERCC-00168" "ERCC-00170" "ERCC-00171"

oGRLgenes.original = oGRLgenes
oGRLgenes = lapply(oGRLgenes, function(gr){
  pass1 = gsub('chr|chrUn_|chr\\d+_', '', as.character(seqnames(gr)))
  pass2 = gsub('_random', '.1', pass1)
  pass3 = gsub('X_|Y_', '', pass2)
  pass4 = gsub('^(\\w+\\d\\d+)$', '\\1.1', pass3)
  gr.new = GRanges(seqnames = pass4,ranges = ranges(gr), strand = strand(gr))
  df = DataFrame(exon_id=gr$exon_id, exon_name=gr$exon_name)
  mcols(gr.new) = df
  return(gr.new)
})

oGRLgenes = GRangesList(oGRLgenes)

# table(seqnames(seqinfo(bf)) %in% seqlevels(oGRLgenes))
# load the ERCC gtf with spikein information
oGRercc = import('~/Data/MetaData/ERCC92.gtf')
strand(oGRercc) = '*'
# reformat metadata column
f = oGRercc$gene_id
df = DataFrame(exon_id=oGRercc$transcript_id, exon_name=NA)
mcols(oGRercc) = df
oGRLercc = split(oGRercc, f)

oGRLgenes = append(oGRLgenes, oGRLercc)


## create the bamfile list from database
## connect to mysql database to get sample information
library('RMySQL')

##### connect to mysql database to get samples
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)
dbListFields(db, 'File')
# get the query
g_did
q = paste0('select Sample.id as sid, Sample.group1, Sample.group2, Sample.title, File.* from Sample, File
           where (Sample.idData = 28) AND (File.idSample = Sample.id AND File.type like "%duplicates removed%")')
dfSample = dbGetQuery(db, q)
nrow(dfSample)
dfSample
# close connection after getting data
dbDisconnect(db)
# remove any whitespace from the names
dfSample$name = gsub(" ", "", dfSample$name, fixed = T)
dfSample$title = gsub(" ", "", dfSample$title, fixed = T)
dfSample$group1 = gsub(" ", "", dfSample$group1, fixed = T)
# create a new file path variable 
dfSample$fp = paste0(dfSample$name)

#### set working directory to appropriate location with bam files
setwd('dataExternal/bamRD/')
csFiles = list.files('.', pattern = '*.bam$', recursive = T)
# check if these files match the file names in database
table(dfSample$fp %in% csFiles)

## order the data frame by sample
dfSample = dfSample[order(dfSample$title),]
# split the files by titles to reduce memory usage
lFiles = split(dfSample$fp, dfSample$title)

# for each of these bam file lists do the counting
lCounts = lapply(lFiles, function(bfl){
  ## create a bamfiles list object
  oBamFiles = BamFileList(bfl, index=paste0(bfl, '.bai'))
  return(assays(summarizeOverlaps(oGRLgenes, oBamFiles, ignore.strand = F, singleEnd=F))$counts)
})


## save the summarized experiment object
# setwd(gcswd)
# n = make.names(paste('lCounts rd id 28 geraldine rds'))
# n2 = paste0('~/Data/MetaData/', n)
# save(lCounts, file=n2)
# 
# ## comment out after first time execution
# library('RMySQL')
# db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
# dbListTables(db)
# dbListFields(db, 'MetaFile')
# df = data.frame(idData=g_did, name=n, type='rds', location='~/Data/MetaData/',
#                 comment='list of Count matrix from geraldine ganoids mouse sequencing run with quality 10 and duplicates removed')
# dbWriteTable(db, name = 'MetaFile', value=df, append=T, row.names=F)
# dbDisconnect(db)
