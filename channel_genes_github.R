# GENES for channel receptors
rm(list=ls())

## INSERT YOUR PATH FILES HERE ##
# Start with an log fold change S4 object as output from DESeq2
path_dir <- '/.../.../DATA/'                     #USER INPUT
DEObject_name  <- 'DESeq2S4_Object'              #USER INPUT
#################################


library(biomaRt)
library(dplyr)
library(hash)


df <- load(file=paste(path_dir, DEObject_name,sep=''))

get_channel_genes <- function(X){
  # X is results dataframe (df)
  # The values of the dict keys are the starting letters for these gene types
  # Channels were shamelessly found on Wikipedia 
  channels <- hash()
  channels[['Potassium']] <- c('KC')
  channels[['Sodium']]    <- c('SCN')
  channels[['Chloride']]  <- c('CLC', 'CLIC', 'BSND', 'CFTR')
  channels[['Calcium']]   <- c('CACN')
  channels[['Glutm']]     <- c('GRM')   # metabotropic glutamate
  channels[['GABAA']]     <- c('GABR')
  channels[['Serotonin']] <- c('HTR')
  channels[['ACh']]       <- c('CHRN')
  channels[['Zinc']]      <- c('ZACN')
  channels[['Glycine']]   <- c('GLR')
  channels[['AMPA']]      <- c('GRIA')
  channels[['KAINATE']]   <- c('GRIK')
  channels[['NMDA']]      <- c('GRIN')
  channels[['Orphan']]    <- c('GRID') 
  channels[['P2X']]       <- c('P2RX')
  channels[['TRP']]       <- c('TRP')  # Transient Receptor
  channels[['TRPV1_Int']] <- c('CALM1','SNAPAP','SYT9','CBD','AEA','NPR1','PKG') #TRVP1 interactions
  
  df_channels <- X[FALSE,]
  grepper <- function(grepee){
    # grepee is the search term
    grepee <- paste('^',grepee,sep='')
    return(X[grep(grepee, rownames(X)),])
  }
  
  for (k in keys(channels)){
    for (v in channels[[k]]){
      grepped <- grepper(v)
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

channels <- get_channel_genes(#Insert name of your dataframe here)    #USER INPUT
GOI_df <- query_ensembl_functions(channels)
write.table(GOI_df, file=paste(path_dir, DEObject_name, 'channels.csv', sep=''), col.names = NA)

