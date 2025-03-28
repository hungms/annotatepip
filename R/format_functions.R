#' convert_list_to_df
#'
#' Store each vector from a list into a dataframe columns
#' @param list A list of vectors
#' @return A dataframe
#' @export
convert_list_to_df <- function(list){
    max_length <- max(sapply(list, length))
    padded.list <- lapply(list, function(v) {
        c(v, rep("", max_length - length(v)))})
    df <- as.data.frame(padded.list, stringsAsFactors = FALSE)
    return(df)}

#' read_gmt
#'
#' Read GMT file from path
#' @param gmt Path to GMT file
#' @return A dataframe
#' @export
read_gmt <- function(gmt){
    stopifnot(file.exists(gmt))
    gmt.df <- read.table(gmt, sep = "\t", header = F, row.names = 1)
    gmt.df <- gmt.df[, -c(1)] %>% t(.) %>% as.data.frame(.)
    return(gmt.df)}

#' read_gct
#'
#' Read GCT file from path
#' @param gmt Path to GCT file
#' @return A dataframe
#' @export
read_gct <- function(gct){
    stopifnot(file.exists(gct))
    gct.df <- read.table(gct, sep = "\t", header = T, row.names = 1, skip = 2)
    gct.df <- gct.df[, -c(1)]
    return(gct.df)}

#' write_gmt
#'
#' Save dataframe as GMT file
#' @param df Dataframe
#' @param save_name File name
#' @param save_dir File path
#' @return Dataframe in GMT format
#' @export
write_gmt <- function(df, save_name = "", save_dir = getwd()){
    gmt.df <- df %>% t(.) %>% as.data.frame(.) %>% rownames_to_column("Name")
    gmt.df$Description <- gmt.df$Name
    gmt.df <- gmt.df[,c(ncol(gmt.df), 1:(ncol(gmt.df)-1))]
    write.table(gmt.df, file = paste0(save_dir, "/", save_name), col.names = F, row.names = F, quote = F, sep = "\t")}

#' write_gct
#'
#' Save dataframe as GCT file
#' @param df Dataframe
#' @param save_name File name
#' @param save_dir File path
#' @return Dataframe in GCT format
#' @export
write_gct <- function(df, save_name = "", save_dir = getwd()){
    gct.df <- matrix(ncol = ncol(df) + 2, nrow = nrow(df) + 3)
    gct.df[1, 1] <- "#1.3"
    gct.df[2, 1:2] <- dim(df)
    gct.df[3, 1:2] <- c("Name", "Description")
    gct.df[3, 3:ncol(gct.df)] <- colnames(df)
    gct.df[4:(nrow(gct.df)), 1] <- rownames(df)
    gct.df[4:(nrow(gct.df)), 2] <- rownames(df)
    gct.df[4:nrow(gct.df), 3:ncol(gct.df)] <- as.matrix(df)
    gct.df <- as.data.frame(gct.df)
    write.table(gct.df, file = paste0(save_dir, "/", save_name), col.names = F, row.names = F, quote = F, sep = "\t")}

#' add_gene_annotations
#'
#' Save dataframe as GCT file
#' @param df Dataframe
#' @param gene_column Column contaning gene symbols
#' @param org Organism, either "human" or "mouse". Default is "human".
#' @param release Ensembl release. Default is "105".
#' @return Dataframe with functional annotation of genes
#' @export
run_annotation <- function(df, gene_column = "gene", org = "human", release = "105"){
    stopifnot(org %in% c("human", "mouse"))
    stopifnot(gene_column %in% colnames(df))
    stopifnot(is.data.frame(df))
    if(org == "mouse"){
        stopifnot(any(str_detect(df[[gene_column]], "[a-z]")) & org == "mouse")}
    else{
        stopifnot(!any(str_detect(df[[gene_column]], "[a-z]")) & org == "human")}

    files <- list.files(system.file("extdata", package = "annotatepip"), full.names = T)
    files <- files[str_detect(files, paste0("release-", release, "_", org, "_omnipath_db.tsv"))]
    annot <- read.table(files, header = T, sep = "\t")
    df <- merge(df, annot, by.x = gene_column, by.y = paste0(org, "_gene_symbol"), all.x = T)
    return(df)
}