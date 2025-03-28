#' get_xy_genes
#'
#' Get gene symbols from the XY chromosome.
#' @param org Organism either "human" or "mouse". Default is "human".
#' @return A vector of XY genes.
#' @export
get_xy_genes <- function(org = "human"){
  stopifnot(org %in% c("mouse", "human"))
  if(!"biomart_dict" %in% ls()){
      import_biomart_local()}
  xy.genes <- biomart_dict %>% 
      filter(!!sym(paste0(org, "_chromosome")) %in% c("X", "Y")) %>%
      dplyr::select(!!sym(paste0(org, "_gene_symbol"))) %>%
      .[[1]] %>%
      unique(.)
    return(xy.genes)
}

#' get_mt_genes
#'
#' Get gene symbols from mitochondrial chromosome.
#' @param org Organism either "human" or "mouse". Default is "human".
#' @return A vector of MT genes.
#' @export
get_mt_genes <- function(org = "human"){
    stopifnot(org %in% c("mouse", "human"))
    if(!"biomart_dict" %in% ls()){
        import_biomart_local()}
    mt.genes <- biomart_dict %>% 
        filter(!!sym(paste0(org, "_chromosome")) %in% c("MT")) %>%
        dplyr::select(!!sym(paste0(org, "_gene_symbol"))) %>%
        .[[1]] %>%
        unique(.)
    return(mt.genes)
}

#' get_str_genes
#'
#' Get gene symbols from default strings.
#' @param org Organism either "human" or "mouse". Default is "human".
#' @param str String to search for. Default is "rb".
#' @return A vector of str-matching genes.
#' @export
get_str_genes <- function(org = "human", str = c("rb")){
    stopifnot(all(str %in% c("bcr", "tcr", "mhc", "hb", "rb", "mt")))
    str <- paste0(unlist(mget(paste0(str, ".string"), envir = getNamespace("annotatepip"))), collapse = "|")

    if(!"biomart_dict" %in% ls()){
        import_biomart_local()}
    genes <- biomart_dict %>% 
        filter(str_detect(!!sym(paste0(org, "_gene_symbol")), str)) %>%
        as.list(.) %>%
        unlist(.) %>%
        unique(.)
    return(genes)}

#' BCR string
#'
#' BCR string
#' @export
bcr.string <- "^I[Gg][HKLhkl][VDJCAEMGLvdjcaemgl]|^AC233755"

#' TCR string
#'
#' TCR string
#' @export
tcr.string <- "^T[Rr][ABCDGabcdg][VDJCvdjc]"

#' HB string
#'
#' HB string
#' @export
hb.string <- "^H[B][ABDEGMPQZ]?\\d*$|^H[b][abdegmpqz]?\\d*"

#' MHC string
#'
#' MHC string
#' @export
mhc.string <- "^HLA-|^H2-"

#' RB string
#'
#' RB string
#' @export
rb.string <- "^R[Pp][SsLl]"

#' MT string
#'
#' MT string
#' @export
mt.string <- "^[Mm][Tt]-"