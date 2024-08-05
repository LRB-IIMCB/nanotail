#' Read Single Nanopolish polyA preditions from file
#'
#' This is the basic function used to import output from \code{nanopolish polya} to R
#'
#' @param polya_path path to nanopolish output file
#' @param sample_name sample name (optional), provided as a string.
#' If specified will be included as an additional column sample_name.
#' @param gencode are contig names GENCODE-compliant.
#' Can get transcript names and ensembl_transcript IDs if reads were mapped for example to Gencode reference transcriptome
#'
#' @seealso \link{read_polya_multiple}
#'
#' @export
#'
#' @return a [tibble][tibble::tibble-package] with polya predictions
#'
read_polya_single <- function(polya_path, gencode = TRUE, sample_name = NA, dorado = FALSE) {
  # required asserts
  
  #check if parameters are provided
  if (missing(polya_path)) {
    stop("The path to polyA predictions (argument polya_path) is missing",
         call. = FALSE)
  }
  assertthat::assert_that(assertive::is_a_non_missing_nor_empty_string(polya_path),msg = "Empty string provided as an input. Please provide a polya_path as a string")
  assertthat::assert_that(assertive::is_existing_file(polya_path),msg=paste("File ",polya_path," not exists",sep=""))
  assertthat::assert_that(assertive::is_non_empty_file(polya_path),msg=paste("File ",polya_path," is empty",sep=""))
  assertthat::assert_that(assertive::is_a_bool(gencode),msg="Please provide TRUE/FALSE values for gencode parameter")
  
  message(paste0("Loading data from ",polya_path))
  
  #integer64 set to "numeric" to avoid inconsistences when called from read_polya_multiple
  
  file_header <- read.table(polya_path,nrows=1)
  if (sum(c("pt","reference","ref_start") %in% file_header)==3) {
    message("Seems like output from dorado. ")
    dorado = TRUE
  }
  
  polya_data <- data.table::fread(polya_path, integer64 = "numeric", data.table = F,header=TRUE,stringsAsFactors = FALSE,check.names = TRUE,showProgress = FALSE) %>% dplyr::as_tibble()
  if (!dorado) {
    polya_data <- polya_data %>% dplyr::mutate(polya_length = round(polya_length),dwell_time=transcript_start-polya_start)
  }
  else {
    polya_data <- polya_data %>% dplyr::mutate(polya_length = round(pt),dwell_time=NA) %>% dplyr::rename(contig=reference)
  }
  # change first column name
  colnames(polya_data)[1] <- "read_id"
  # transcript names, if mapping to gencode transcriptome
  if (gencode == TRUE) {
    transcript_names <- gsub(".*?\\|.*?\\|.*?\\|.*?\\|.*?\\|(.*?)\\|.*", "\\1", polya_data$contig)
    polya_data$transcript <- transcript_names
    ensembl_transcript_ids <- gsub("^(.*?)\\|.*\\|.*", "\\1", polya_data$contig)
    ensembl_transcript_ids_short <- gsub("(.*)\\..*", "\\1", ensembl_transcript_ids) # without version number
    polya_data$ensembl_transcript_id_full <- ensembl_transcript_ids
    polya_data$ensembl_transcript_id_short <- ensembl_transcript_ids_short
  }
  else {
    polya_data <- polya_data %>% dplyr::rename(transcript = contig)
  }
  
  if(!is.na(sample_name)) {
    # set sample_name (if was set)
    if (! "sample_name" %in% colnames(polya_data)) {
      warning("sample_name was provided in the input file. Overwriting with the provided one")
    }
    polya_data$sample_name = sample_name
    polya_data$sample_name <- as.factor(polya_data$sample_name)
  }
  
  return(polya_data)
}



#' Read Single Nanopolish polyA preditions from file
#'
#' This is the basic function used to import output from poly(A) prediction programs to R
#'
#' @param polya_path path to nanopolish output file
#' @param sample_name sample name (optional), provided as a string.
#' If specified will be included as an additional column sample_name.
#' @param gencode are contig names GENCODE-compliant.
#' Can get transcript names and ensembl_transcript IDs if reads were mapped for example to Gencode reference transcriptome
#' @param input_type type of input file. Can be one of: \itemize{
#' \item nanopolish - nanopolish output
#' \item dorado - dorado output
#' \item tailfindr - tailfindr output
#' \item auto - automatically detect input type
#' }
#' 
#' @seealso \link{read_polya_multiple}
#'
#' @export
#'
#' @return a [tibble][tibble::tibble-package] with polya predictions
#'
read_polya_single2 <- function(polya_path, gencode = TRUE, sample_name = NA, input_type = "auto") {
    # required asserts

    #check if parameters are provided
    if (missing(polya_path)) {
      stop("The path to polyA predictions (argument polya_path) is missing",
           call. = FALSE)
    }
    

    # checkmate - input_type in values: nanopolish, dorado, tailfindr, auto
    checkmate::assert_string(input_type)
    checkmate::assert_choice(input_type, c("nanopolish", "dorado", "tailfindr","auto","dorado_old"))

    checkmate::assert_string(polya_path)
    checkmate::assert_file_exists(polya_path)
    checkmate::assert_logical(gencode)

    

    # assert the file in polya_path is not empty:
    if (file.size(polya_path) == 0) {
      stop("File ",polya_path," is empty",call. = FALSE)
    }
  
  
  
    message(paste0("Loading data from ",polya_path))

   
    
    
    file_header <- data.table::fread(polya_path,nrows=0,header=T,data.table=F) # read first line to check colnames to determine input type
    message(file_header)
    # auto-detect input type
    
    if ((sum(c("reference","ref_start","pt") %in% colnames(file_header))==3)) {
        guessed_input_type <- "dorado_old" # for compatibility with previous versions of python script getting poly(A) from bams 
        message("Dorado output detected (older naming convention with pt insted of polya_length column)")
    }
    else if ((sum(c("reference","ref_start") %in% colnames(file_header))==2)) {
      guessed_input_type <- "dorado"
      message("Dorado output detected")
    }
    else if (sum(c("read_type","tail_is_valid") %in% colnames(file_header))==2) {
      guessed_input_type <- "tailfindr"
      message("Tailfindr output detected")
    }
    else if (sum(c("read_rate","qc_tag") %in% colnames(file_header))==2) {
      guessed_input_type <- "nanopolish"
      message("Nanopolish output detected")
    }
    else {
        stop("Unknown type of input file",call. = FALSE)
    }
    
    if (input_type!="auto" & (guessed_input_type != input_type)) {
        warning("Input type does not match the detected one. Detected: ",guessed_input_type,". Specified: ",input_type)
    }
    
    if (input_type=="auto") {
        input_type <- guessed_input_type
    }
    
    if (input_type=="nanopolish") {
      columns_to_keep = c("readname","contig","position","polya_length","qc_tag")
    }
    else if (input_type=="dorado") {
      columns_to_keep = c("read_id","reference","ref_start","polya_length","mapq")
    }
    else if (input_type=="dorado_old") {
      columns_to_keep = c("read_id","reference","ref_start","pt","mapq")
    }
    else if (input_type=="tailfindr") {
      columns_to_keep = c("read_id","read_type","tail_length","tail_is_valid")
    }

    
    
    #integer64 set to "numeric" to avoid inconsistences when called from read_polya_multiple
    polya_data <- data.table::fread(polya_path, integer64 = "numeric", data.table = F,header=TRUE,stringsAsFactors = FALSE,check.names = TRUE,showProgress = FALSE, select = columns_to_keep) %>% dplyr::as_tibble()
    
    # create consistent output, with the same column names and order for each input type
    # colnames: c("readname","reference","ref_start","polya_length","qc_tag","mapq","read_type")
   
    
    if (input_type=='dorado' | input_type=='dorado_old') {
      polya_data$qc_tag <- "PASS"
      polya_data$read_type <- NA
      # make proper order of columns
      polya_data <- polya_data[,c(1,2,3,4,6,5,7)]
      colnames(polya_data) <- c("readname","reference","ref_start","polya_length","qc_tag","mapq","read_type")
    }    
    else if (input_type=='nanopolish') {
      polya_data$mapq <- NA
      polya_data$read_type <- NA
      colnames(polya_data) <- c("readname","reference","ref_start","polya_length","qc_tag","mapq","read_type")
    }
    else if (input_type=='tailfindr') {
      polya_data$reference <- NA
      polya_data$ref_start <- NA
      polya_data$mapq <- NA
      # make reference second column and ref_start third column in polya_data using base R
      polya_data <- polya_data[,c(1,5,6,3,4,7,2)]
      colnames(polya_data) <- c("readname","reference","ref_start","polya_length","qc_tag","mapq","read_type")
      polya_data$qc_tag <- ifelse(polya_data$qc_tag==TRUE,"PASS","FAIL")
    }
    
    # change first column name
    colnames(polya_data)[1] <- "read_id"
    # transcript names, if mapping to gencode transcriptome
    # skip if tailfindr output - as it lacks information on the reference
    if (gencode == TRUE & input_type!="tailfindr") {
        transcript_names <- gsub(".*?\\|.*?\\|.*?\\|.*?\\|.*?\\|(.*?)\\|.*", "\\1", polya_data$reference)
        polya_data$transcript <- transcript_names
        ensembl_transcript_ids <- gsub("^(.*?)\\|.*\\|.*", "\\1", polya_data$reference)
        ensembl_transcript_ids_short <- gsub("(.*)\\..*", "\\1", ensembl_transcript_ids) # without version number
        polya_data$ensembl_transcript_id_version <- ensembl_transcript_ids
        polya_data$ensembl_transcript_id <- ensembl_transcript_ids_short
    }
    else {
      # if not gencode use contig (mapped reference) as transcript name
      polya_data$transcript <- polya_data$reference
    }

    if(!is.na(sample_name)) {
      # set sample_name (if was set)
      if ("sample_name" %in% colnames(polya_data)) {
        warning("sample_name was provided in the input file. Overwriting with the provided one")
      }
      polya_data$sample_name = sample_name
      polya_data$sample_name <- as.factor(polya_data$sample_name)
    }

    return(polya_data)
}


# TODO Benchmark na lapply/for loop

#' Reads multiple polyA predictions at once
#'
#' This function can be used to load any number of files with polyA predictions with single invocation,
#' allowing for metadata specification.
#'
#'
#' @param samples_table data.frame or tibble containing samples metadata and paths to files.
#' Should have at least two columns: \itemize{
#' \item polya_path - containing path to the polya predictions file
#' \item sample_name - unique name of the sample
#' }
#' Additional columns can provide metadata which will be included in the final table
#' @param ... - additional parameters to pass to read_polya_single(), like gencode=(TRUE/FALSE)
#'
#' @return a [tibble][tibble::tibble-package] containing polyA predictions for all specified samples, with metadata provided in samples_table
#' stored as separate columns
#'
#' @seealso \link{read_polya_single}
#'
#' @export
#'
read_polya_multiple <- function(samples_table,...) {

  if (missing(samples_table)) {
    stop("Samples table argument is missing",
         call. = FALSE)
  }

  assertthat::assert_that(assertive::has_rows(samples_table),msg = "Empty data frame provided as an input (samples_table). Please provide samples_table describing data to load")
  assertthat::assert_that("polya_path" %in% colnames(samples_table),msg = "Samples table should contain at least polya_path and sample_name columns")
  assertthat::assert_that("sample_name" %in% colnames(samples_table),msg = "Samples table should contain at least polya_path and sample_name columns")

  samples_data <- samples_table %>% dplyr::as.tbl() %>% dplyr::mutate_if(is.character,as.factor) %>% dplyr::mutate(polya_path = as.character(polya_path)) %>% dplyr::group_by(sample_name) %>% dplyr::mutate(polya_contents=purrr::map(polya_path, function(x) read_polya_single(x))) %>% dplyr::ungroup() %>% dplyr::select(-polya_path)
  polya_data <- tidyr::unnest(samples_data)

  return(polya_data)
}


#' Removes reads which failed during Nanopolish polya processing
#'
#' Convenient function to quickly remove all reads failing during nanopolish polya processing
#'
#' @param polya_data output table from \link{read_polya_single} or \link{read_polya_multiple}
#'
#' @return a [tibble][tibble::tibble-package] with only reads having qc_tag=='PASS'
#'
#' @export
#'
#' @seealso \link{read_polya_single}, \link{read_polya_multiple}
#'
remove_failed_reads <- function(polya_data) {

  if (missing(polya_data)) {
    stop("Please provide data.frame with polyA predictions as an input.",
         call. = FALSE)
  }

  assertthat::assert_that(assertive::has_rows(polya_data),msg = "Empty data frame provided as an input (polya_data). Please provide valid input")

  filtered_polya_data <- polya_data %>% dplyr::filter(qc_tag=='PASS')
  return(filtered_polya_data)
}


