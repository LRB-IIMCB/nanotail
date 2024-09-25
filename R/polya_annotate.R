#' Annotate polyA predictions using annotables
#'
#' @param polya_data polya data table to annotate
#' @param genome valid genome from annotables to use for annotation
#'
#' @return a \link[tibble]{tibble}
#' @export
#'
annotate_with_annotables <- function(polya_data,genome) {

  if ( !requireNamespace('annotables',quietly = TRUE) ) {
    stop("NanoTail requires 'annotables'. Please install it using
         install.packages('devtools')
        devtools::install_github('stephenturner/annotables')")
  }
  require(annotables)


  if (missing(polya_data)) {
    stop("Please provide data.frame with polyA predictions as an input.",
         call. = FALSE)
  }

  if (missing(genome)) {
    stop("Please provide valid genome from annotables package to use for annotation.",
         call. = FALSE)
  }

  # assertthat::assert_that(assertive::has_rows(polya_data),msg = "Empty data frame provided as an input (polya_data). Please provide valid input")

  tx_to_gene_table = paste0(genome,"_tx2gene")
  # Annotatio join, with last step to deduplicate annotations - remove duplicates which occurs due to e.g multiple entrez ids for each transcript
  polya_data_annotated <-  polya_data %>% dplyr::left_join( eval(as.symbol(tx_to_gene_table)),by=c("ensembl_transcript_id_short"  = "enstxp")) %>% dplyr::left_join(eval(as.name(genome)) %>% dplyr::group_by(ensgene) %>% dplyr::slice(1) %>% dplyr::ungroup())

  #Explictly convert selected columns to factors (required for proper visualization)
  polya_data_annotated$biotype <- as.factor(polya_data_annotated$biotype)
  polya_data_annotated$strand <- as.factor(polya_data_annotated$strand)
  polya_data_annotated$chr <- as.factor(polya_data_annotated$chr)




  return(polya_data_annotated)
  }


#' Title
#'
#' @param polya_data polya data table to annotate
#' @param attributes_to_get what annotations should be retrieved. Default = c('external_gene_name','description','transcript_biotype')
#' @param filters which column should be matched in the target mart
#' @param mart_to_use mart object created with \link[biomaRt]{useMart} or \link[biomaRt]{useEnsembl}
#'
#' @return a \link[tibble]{tibble}
#' @export
annotate_with_biomart <- function(polya_data,attributes_to_get=c('ensembl_transcript_id','external_gene_name','description','transcript_biotype'),filters='ensembl_transcript_id',mart_to_use=NA) {

  if (missing(polya_data)) {
    stop("Please provide data.frame with polyA predictions as an input.",
         call. = FALSE)
  }

  if (missing(mart_to_use)) {
    stop("Please provide valid mart object",
         call. = FALSE)
  }

  # assertthat::assert_that(assertive::has_rows(polya_data),msg = "Empty data frame provided as an input (polya_data). Please provide valid input")
  assertthat::assert_that(class(mart_to_use)=="Mart",msg="Please provide valid mart object")
  assertthat::assert_that(length(attributes)>0,msg="please provide attributes")

  ensembl_ids = unique(polya_data$ensembl_transcript_id_short)
  ensembl_ids <- ensembl_ids[!is.na(ensembl_ids)]

  polya_data <- polya_data %>% dplyr::rename(ensembl_transcript_id = ensembl_transcript_id_short)

  # if using biomaRt version older than from Bioconductor 3.9, it cannot process more than 500 values at once
  if (packageVersion("biomaRt")<"2.40.0") {
    number_of_items <- length(ensembl_ids)
    annotation_data=data.frame()
    for (z in seq(1,number_of_items,by = 500)) {

      annotation_data_temp<-getBM(attributes=attributes_to_get, filters =filters, values = ensembl_ids[z:(z+499)], mart = mart_to_use)
      #print(annotation_data_temp)
      annotation_data<-rbind(annotation_data,annotation_data_temp)
    }
  }
  # since biomaRt 2.40 batch submission is possible
  else {
    annotation_data<-biomaRt::getBM(attributes=attributes_to_get, filters = filters, values = ensembl_ids, mart = mart_to_use)
  }
  polya_data_annotated <-  polya_data %>% dplyr::left_join(annotation_data)

  return(polya_data_annotated)
}


#' Title
#'
#' @param columns_of_annotation which columns to use
#' @param keytype whic keytype to use
#' @param organism whic organism database to use
#' @param polya_data polya data table to annotate
#'
#' @return a \link[tibble]{tibble}
#' @export
annotate_with_org_packages <- function(polya_data,columns_of_annotation=c("GENENAME","SYMBOL"),keytype='ENSEMBLTRANS',organism="mus_musculus") {

  if (missing(polya_data)) {
    stop("Please provide data.frame with polyA predictions as an input.",
         call. = FALSE)
  }


  # currently thos supported
  valid_org_packages = list("homo_sapiens" = "org.Hs.eg.db", "mus_musculus" = "org.Mm.eg.db","rattus_norvegicus" = "org.rn.eg.db","saccharomyces_cerevisiae" = "org.Sc.sgd.db","caenorhabditis_elegans" = "org.Ce.eg.db")

  # assertthat::assert_that(assertive::has_rows(polya_data),msg = "Empty data frame provided as an input (polya_data). Please provide valid input")
  assertthat::assert_that(length(columns_of_annotation)>0,msg="please provide columns of annotation")

  ensembl_ids = unique(polya_data$ensembl_transcript_id_short)
  ensembl_ids <- ensembl_ids[!is.na(ensembl_ids)]



  polya_data <- polya_data %>% dplyr::rename(!! rlang::sym(keytype) := ensembl_transcript_id_short)




  annotation_data<-AnnotationDbi::select(eval(parse(text = valid_org_packages[[organism]])),columns = columns_of_annotation,keytype = keytype,keys = ensembl_ids)

  polya_data_annotated <-  polya_data %>% dplyr::left_join(annotation_data)

  return(polya_data_annotated)
}




#' Parse references provided in the gencode format 
#' 
#' 
#' @param input_vector vector with gencode-formatted references
#' 
#' @return a \link[tibble]{tibble}
#' @export
#' 
#' @examples
#' parse_gencode_headers(c("ENST00000456328|ENSG00000223997|OTTHUMG00000000961|OTTHUMT00000002844.2|OTTHUMT00000002844|DDX11L1|202|processed_transcript|","ENST00000450305|ENSG00000223997|OTTHUMG00000000961|OTTHUMT00000002844.2|OTTHUMT00000002844|DDX11L1|202|processed_transcript|"))
#' 
#' 
parse_gencode_headers <- function(input_vector) {
  
  # check if input_list is provided
  if (missing(input_vector)) {
    stop("Input vector with gencode-formatted references is missing. Please provide a valid argument",
         call. = FALSE)
  }
  
  # check if input_list has at least one element
  if (length(input_vector) < 1) {
    stop("Input vector should have at least one element",
         call. = FALSE)
  }
  
  
  gencode_elements_no = sum(grepl("^ENST.*\\|ENSG.*\\|$",input_vector))
  vector_length = length(input_vector)
  
  
  if (!gencode_elements_no>0) {
    stop("No gencode reference provided",
         call. = FALSE)
  }
  
  # check if input_vector 1st element is formatted as gencode transcriptome reference
  if (!gencode_elements_no==vector_length) {
    warning("Not all of provided vector elements contain gencode-formatted reference. Proceeding anyway")
  }
  
  
  extract_annotations <- function(x) {
    parts <- strsplit(x, "\\|")[[1]]
    ensembl_transcript_id <- ifelse(length(parts) >= 1, parts[1], x)
    ensembl_gene_id <- ifelse(length(parts) >= 3, parts[3], x)
    symbol <- ifelse(length(parts) >= 6, parts[6], x)
    return(c(reference=x,ensembl_transcript_id = ensembl_transcript_id, ensembl_gene_id = ensembl_gene_id, symbol = symbol))
  }
  
  annotations <- data.frame(t(sapply(input_vector, extract_annotations,USE.NAMES = FALSE)))
  
  return(annotations)
  
}