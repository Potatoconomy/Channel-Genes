# GENES for channel receptors
rm(list=ls())

## INSERT YOUR PATH FILES HERE ##
# Start with an log fold change S4 object as output from DESeq2
path_dir <- '/YOUR/DATA/DIR/'         #USER INPUT
DEObject_name  <- 'DESeq2_Object'     #USER INPUT
#################################


library(biomaRt)
library(dplyr)
library(hash)


df <- load(file=paste(path_dir, DEObject_name,sep=''))

get_channel_genes <- function(X){
  # X is results dataframe (df)
  # The values of the dict keys are the starting letters for these gene types
  # Channels were shamelessly found on Wikipedia, lab member advice, GO, and genenames.org 
  channels <- hash()
  channels[['Potassium']]         <- c('KCN')  #https://www.genenames.org/data/genegroup/#!/group/183
  channels[['Sodium']]            <- c('SCN')
  channels[['Chloride']]          <- c('CLC', 'CLIC', 'BSND', 'CFTR')
  channels[['Calcium']]           <- c('CACN')
  channels[['Glutm']]             <- c('GRM')   # metabotropic glutamate
  channels[['GABAA']]             <- c('GABR')
  channels[['Serotonin']]         <- c('HTR')
  channels[['ACh']]               <- c('CHRN')
  channels[['Zinc']]              <- c('ZACN')
  channels[['Glycine']]           <- c('GLR')
  channels[['AMPA']]              <- c('GRIA')
  channels[['KAINATE']]           <- c('GRIK')
  channels[['NMDA']]              <- c('GRIN')
  channels[['Orphan']]            <- c('GRID') 
  channels[['P2X']]               <- c('P2RX')
  channels[['TRP']]               <- c('TRP')  # Transient Receptor
  channels[['TRPV1_Int']]         <- c('CALM1','SNAPAP','SYT9','CBD','AEA','NPR1','PKG') #TRVP1 interactions
  channels[['glutamatesynapse']]  <- c("SYT1","SYN2","VAMP2","BSN","SLC17A6","SLC17A7", "SLC17A8")
  channels[['Presynaptic']]       <- c("SYT1","SYN2","VAMP2","BSN")
  channels[['Postsynaptic']]      <- c("DLG4","NLGN1","NLGN2","HOMER1")
  channels[['vGLUT']]             <- c("SLC17A6","SLC17A7", "SLC17A8")
  channels[['astrosynapse']]      <- c("GLUL","SLC38A3","GRM3","GRM5","SLC1A2","SLC1A3")
  channels[['Mglur']]             <- c("GRM1","GRM2","GRM3","GRM4","GRM5","GRM6","GRM7","GRM8")
  channels[['AIS']]               <- c("ANK3","BCAN","BIN1","CAMK2D","CCK","CNGA3","IQSCHFP","KCNQ2",
                                       "KCNQ3","LRRC7","MAP1A","MAP2","NAV1","NFASC","NRCAM","SCN1A",
                                       "SCN2A","SCN8A","SPTBN4","SPTBN","TRIM46")
  channels[['Calciumbinding']]    <- c("CALB1","CALB2")
  channels[['calciumsensor']]     <- c("EFCAB9","EFHB","EFHD1","SYT1")
  channels[['tempsensor']]        <- c("ADORA1","ANO1","ARRB2","ASIC3","CALCA","CXCL12","CXCR4",
                                       "EPHB1","HTR2A","LXN","MMP24","NR2F6","NTRK1","NTSR1",
                                       "PRDM12","TRPV1")
  channels[['synaptogamins']]     <-c("SYT1","SYT2","SYT3","SYT4","SYT5","SYT6","SYT7","SYT8","SYT9",
                                      "SYT10","SYT11","SYT12","SYT13","SYT14","SYT15","SYT16","SYT17")
  channels[['calciumbindingsynaptogamin']] <-c("SYT1","SYT2","SYT3","SYT5","SYT6","SYT7","SYT9")
  channels[['presynapticcell']]   <-c("SLC17A6","SLC17A7", "SLC17A8","PRKACA","PRKACB","PRKACG",
                                      "ADCY1","ADCY2","ADCY3","ADCY5","ADCY6","ADCY7","ADCY8",
                                      "ADCY9","ADCY4","CACNA1A","KCNJ3")
  channels[['postsynapticcell']]  <-c("GRIK","GRIA","GRIN","TRPC1","DLG4","DLGAP1","SHANK1",
                                      "SHANK2","SHANK3","HOMER1","HOMER2","HOMER3","PPP3CA",
                                      "PPP3CB","PPP3CC","PPP3R1","PPP3R2","CACNA1C","CACNA1D")
  
  X$type <- NA
  df_channels <- X[FALSE,]
  grepper <- function(grepee){
    # grepee is the search term
    grepee <- paste('^',grepee,sep='')
    return(X[grep(grepee, rownames(X)),])
  }
  
  for (k in keys(channels)){
    for (v in channels[[k]]){
      grepped <- grepper(v)
      if (nrow(grepped) != 0){
        grepped$type <- k
      }
      print(grepped)
      df_channels <- data.frame(rbind(as.matrix(df_channels), as.matrix(grepped)))
    }
  }
  
  remove_bad_greps <- df_channels[grep('^GLRX', rownames(df_channels)),]
  df_channels <- df_channels[!rownames(df_channels) %in% rownames(remove_bad_greps),]
  
  # To remove future bad greps, just double the last 2 lines but modified
  
  return(df_channels)
}


query_ensembl_functions <- function(df) {
  # input a dataframe of genes of interest   
  #IMPORTANT: DF ROWNAMES MUST BE IN HGNC SYMBOL, OTHERWISE CHANGE ATTRIBUTES
  # Returns basic gene descriptions
  ensembl <- useMart("ensembl")
  ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl) # can change organism here
  
  ensembl <- useMart("ensembl")
  ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
  query <- getBM(attributes=c('hgnc_symbol','description'),
                 mart=ensembl,
                 values=rownames(df),
                 filters='hgnc_symbol',
                 useCache = FALSE)
  
  a <- distinct(as.data.frame(query), hgnc_symbol, .keep_all = TRUE)
  
  idx <- match(rownames(df), query$hgnc_symbol)
  df$description <- query$description[idx]
  
  return(df)
}

channels <- get_channel_genes(#YOUR RESULTS HERE)     # USER INPUT HERE. PUT NAME OF DESEQ2 OBJECT IN
GOI_df <- query_ensembl_functions(channels)
write.table(GOI_df, file=paste(path_dir, DEObject_name, 'channels.csv', sep=''), col.names = NA)