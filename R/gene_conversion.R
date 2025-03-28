#' get_biomart
#'
#' Retrieves archived ensembl biomarts from web.
#' 
#' @param host URL to retrieve archived ensembl biomarts
#' @return biomart_dict in .GlobalEnv
#' @export
import_biomart <- function(host = 'https://dec2021.archive.ensembl.org') {
   assign("human_biomart", biomaRt::useMart("ensembl", host = host, dataset = "hsapiens_gene_ensembl"), envir = .GlobalEnv)
   assign("mouse_biomart", biomaRt::useMart("ensembl", host = host, dataset = "mmusculus_gene_ensembl"), envir = .GlobalEnv)
   biomart_dict <- getLDS(attributes = c("mgi_symbol", "chromosome_name"), # Mouse attributes
                    filters = "", # No specific filter, retrieve all genes
                    values = NULL, # No specific values, retrieve all genes
                    mart = mouse_biomart, # Mouse BioMart
                    attributesL = c("hgnc_symbol", "chromosome_name"), # Human attributes
                    martL = human_biomart, # Human BioMart
                    uniqueRows = TRUE) # Remove duplicate rows
   biomart_dict <- biomart_dict %>% 
    mutate(
        mouse_gene_symbol = MGI.symbol,
        human_gene_symbol = HGNC.symbol,
        mouse_chromosome = Chromosome.scaffold.name,
        human_chromosome = Chromosome.scaffold.name.1) %>%
    dplyr::select(mouse_gene_symbol, mouse_chromosome, human_gene_symbol, human_chromosome) %>%
    filter(str_detect(mouse_gene_symbol, "[a-z]") | str_detect(human_gene_symbol, "[A-Z]"))
   assign("biomart_dict", biomart_dict, envir = .GlobalEnv)   
   }


#' import_biomart_local
#'
#' Retrieves archived ensembl biomarts from local directory.
#' 
#' @param release Ensembl release number
#' @return biomart_dict in .GlobalEnv
#' @export
import_biomart_local <- function(release = "105") {
    file <- paste0(system.file("extdata", package = "annotatepip"), "/release-", release, ".tsv")
    biomart_dict <- read.table(file, header = TRUE, sep = "\t") %>%
      filter(str_detect(mouse_gene_symbol, "[a-z]") & str_detect(human_gene_symbol, "[A-Z]"))
    assign("biomart_dict", biomart_dict, envir = .GlobalEnv)}

#' convert_mouse_to_human
#' 
#' Convert mouse (MGI) gene symbols to human (HGNC) gene symbols.
#' 
#' @param genes A vector of mouse gene symbols.
#' @param one.to.many Enable many human gene symbols to mapped to same mouse gene
#' @return A vector of unique human gene symbols that are mapped one-to-one. If one.to.many = TRUE, returns a dataframe that mapped one-to-many.
#' @export
convert_mouse_to_human <- function(genes, one.to.many = F){
   if(!"biomart_dict" %in% ls()){
      import_biomart_local()}
   biomart_dict_tmp <- biomart_dict %>% 
      dplyr::select(c(human_gene_symbol, mouse_gene_symbol)) %>% 
      distinct(human_gene_symbol, .keep_all = T)
   if(one.to.many){
      return(biomart_dict_tmp)}
   else{
      human.genes <- biomart_dict_tmp$human_gene_symbol %>% unique(.)
      return(human.genes)}
   }

#' convert_human_to_mouse
#' 
#' Convert human (HGNC) gene symbols to mouse (MGI) gene symbols.
#' 
#' @param genes A vector of human gene symbols.
#' @param one.to.many Enable many mouse gene symbols to mapped to same humangene
#' @return A vector of unique mouse gene symbols that are mapped one-to-one. If one.to.many = TRUE, returns a dataframe that mapped one-to-many.
#' @export
convert_human_to_mouse <- function(genes, one.to.many = F){
   if(!"biomart_dict" %in% ls()){
      import_biomart_local()}
   biomart_dict_tmp <- biomart_dict %>% 
      dplyr::select(c(human_gene_symbol, mouse_gene_symbol)) %>% 
      distinct(mouse_gene_symbol, .keep_all = T)
   if(one.to.many){
      return(biomart_dict_tmp)}
   else{
      mouse.genes <- biomart_dict_tmp$mouse_gene_symbol %>% unique(.)
      return(moue.genes)}
   }

