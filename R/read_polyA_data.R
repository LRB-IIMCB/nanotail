#' Read Single Nanopolish polyA preditions from file
#'
#' This is the basic function used to import output from poly(A) prediction programs to R
#'
#' @param polya_path path to nanopolish output file
#' @param sample_name sample name (optional), provided as a string.
#' If specified will be included as an additional column sample_name.
#' @param input_type type of input file. Can be one of: \itemize{
#' \item nanopolish - nanopolish output
#' \item dorado - dorado output
#' \item tailfindr - tailfindr output
#' \item auto - automatically detect input type
#' }
#' @param metadata additional columns to keep in the output table
#' @param additional_columns_to_keep additional columns to keep in the output table
#' 
#' @seealso \link{read_polya_multiple}
#'
#' @export
#'
#' @return a [tibble][tibble::tibble-package] with polya predictions
#'
read_polya_single <- function(polya_path, gencode = TRUE, verbose=TRUE, sample_name = NA, input_type = "auto", metadata,additional_columns_to_keep) {
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

    # assert the file in polya_path is not empty:
    if (file.size(polya_path) == 0) {
      stop("File ",polya_path," is empty",call. = FALSE)
    }
  
    if (!missing(metadata)) {
      message("metadata provided as an argument")
      
    }
  
    if (verbose) {
      message(paste0("Loading data from ",polya_path))
    }
    
   
    
    
    file_header <- data.table::fread(polya_path,nrows=0,header=T,data.table=F) # read first line to check colnames to determine input type
    #message(file_header)
    # auto-detect input type
    
    if ((sum(c("reference","ref_start","pt") %in% colnames(file_header))==3)) {
        guessed_input_type <- "dorado_old" # for compatibility with previous versions of python script getting poly(A) from bams 
        if (verbose) {
          message("Dorado output detected (older naming convention with pt insted of polya_length column)")
        }
    }
    else if ((sum(c("reference","ref_start") %in% colnames(file_header))==2)) {
      guessed_input_type <- "dorado"
      if (verbose) {
        message("Dorado output detected")
      }
    }
    else if (sum(c("read_type","tail_is_valid") %in% colnames(file_header))==2) {
      guessed_input_type <- "tailfindr"
      if (verbose) {
        message("Tailfindr output detected")
      }
    }
    else if (sum(c("read_rate","qc_tag") %in% colnames(file_header))==2) {
      guessed_input_type <- "nanopolish"
      if (verbose) {
        message("Nanopolish output detected")
      }
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

    if (!missing(additional_columns_to_keep)) {
      columns_to_keep <- c(columns_to_keep,additional_columns_to_keep)
      no_additional_columns <- length(additional_columns_to_keep)
      if (verbose) {
        message(paste0("additional_columns: ",no_additional_columns))
        message(columns_to_keep)
      }
    }
    
    #integer64 set to "numeric" to avoid inconsistences when called from read_polya_multiple
    polya_data <- data.table::fread(polya_path, integer64 = "numeric", data.table = F,header=TRUE,stringsAsFactors = FALSE,check.names = TRUE,showProgress = FALSE, select = columns_to_keep) %>% dplyr::as_tibble()
    
    # create consistent output, with the same column names and order for each input type
    # colnames: c("readname","reference","ref_start","polya_length","qc_tag","mapq","read_type")
   
    
    if (input_type=='dorado' | input_type=='dorado_old') {
      polya_data$qc_tag <- "PASS"
      polya_data$read_type <- NA
      # make proper order of columns
      if (!missing(additional_columns_to_keep)) {
        polya_data <- polya_data[,c(1,2,3,4,(no_additional_columns+6),5,(no_additional_columns+7),seq(6,(no_additional_columns+5)))]
      }
      else {
        polya_data <- polya_data[,c(1,2,3,4,6,5,7)]
      }
    }    
    else if (input_type=='nanopolish') {
      polya_data$mapq <- NA
      polya_data$read_type <- NA
      if (!missing(additional_columns_to_keep)) {
        polya_data <- polya_data[,c(1,2,3,4,5,(no_additional_columns+6),(no_additional_columns+7),seq(6,(no_additional_columns+5)))]
      }
      else {
        polya_data <- polya_data[,c(1,2,3,4,6,5,7)]
      }
      #colnames(polya_data) <- c("readname","reference","ref_start","polya_length","qc_tag","mapq","read_type")
    }
    else if (input_type=='tailfindr') {
      polya_data$reference <- NA
      polya_data$ref_start <- NA
      polya_data$mapq <- NA
      # make reference second column and ref_start third column in polya_data using base R
      if (!missing(additional_columns_to_keep)) {
        polya_data <- polya_data[,c(1,5,(no_additional_columns+6),3,4,(no_additional_columns+7),2,seq(5,(no_additional_columns+4)))]
      }
      else {
        polya_data <- polya_data[,c(1,5,6,3,4,7,2)]
      }
        #colnames(polya_data) <- c("readname","reference","ref_start","polya_length","qc_tag","mapq","read_type")
      polya_data$qc_tag <- ifelse(polya_data$qc_tag==TRUE,"PASS","FAIL")
    }
    
    if (!missing(additional_columns_to_keep)) {
      colnames(polya_data) <- c("read_id","reference","ref_start","polya_length","qc_tag","mapq","read_type",additional_columns_to_keep)
    }
    else {
      colnames(polya_data) <- c("read_id","reference","ref_start","polya_length","qc_tag","mapq","read_type")
    }
    # transcript names, if mapping to gencode transcriptome
   

    # if(!is.na(sample_name)) {
    #   # set sample_name (if was set)
    #   if ("sample_name" %in% colnames(polya_data)) {
    #     warning("sample_name was provided in the input file. Overwriting with the provided one")
    #   }
    #   polya_data$sample_name = sample_name
    #   polya_data$sample_name <- as.factor(polya_data$sample_name)
    # }

    return(polya_data)
}




# with the input table, take the path column and use it to read the content of files paths provided in ths column
# return the content of the files as a list
# the name of the list elements should be the same as the sample_name column
# each element should contain the data element with the content read from the file
# another element of each list element should be named meta, and contain the metadata of the sample, taken from all remaining columns of the input table
# use base R function, dplyr or purrr usage is prohibited

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
#'@examples
#' \dontrun{
#' 
#' read_polya_multiple(example_sample_table)
#' 
#' }
#' 
read_polya_multiple <- function(input_table,verbose=TRUE,process_references=TRUE,...) {
  # check if input_table is provided
  if (missing(input_table)) {
    stop("Table is missing. Please provide a valid table argument",
         call. = FALSE)
  }
  # check if input_table is a data.frame
  checkmate::assert_data_frame(input_table)
  # if (!is.data.frame(input_table)) {
  #   stop("Table should be provided as a data.frame",
  #        call. = FALSE)
  # }
  # 
  
  # check if input_table has at least two columns
  if (ncol(input_table) < 2) {
    stop("Table should have at least two columns",
         call. = FALSE)
  }
  
  # check if input_table has a column named polya_path
  if (!("polya_path" %in% colnames(input_table))) {
    stop("Table should have a column named polya_path",
         call. = FALSE)
  }
  
  # check if input_table has a column named sample_name
  if (!("sample_name" %in% colnames(input_table))) {
    stop("Table should have a column named sample_name",
         call. = FALSE)
  }
  
  # check if input_table has a column named sample_id
  if (!("sample_id" %in% colnames(input_table))) {
    stop("Table should have a column named sample_id",
         call. = FALSE)
  }
  
  # check if input_table has a column named polya_path which is character
  if (!all(sapply(input_table$polya_path, is.character))) {
    stop("Column polya_path should be character",
         call. = FALSE)
  }
  
  # check if input_table has a column named sample_name which is character
  if (!all(sapply(input_table$sample_name, is.character))) {
    stop("Column sample_name should be character",
         call. = FALSE)
  }
  # check if input_table has a column named sample_id which is character
  if (!all(sapply(input_table$sample_id, is.character))) {
    stop("Column sample_id should be character",
         call. = FALSE)
  }
  
  # check if input_table has a column named polya_path which is not empty
  if (!all(sapply(input_table$polya_path, function(x) nchar(x) > 0))) {
    stop("Column polya_path should not be empty",
         call. = FALSE)
  }
  
  # check if input_table has a column named sample_name which is not empty
  if (!all(sapply(input_table$sample_name, function(x) nchar(x) > 0))) {
    stop("Column sample_name should not be empty",
         call. = FALSE)
  }
  
  # check if input_table has a column named sample_id which is not empty
  if (!all(sapply(input_table$sample_id, function(x) nchar(x) > 0))) {
    stop("Column sample_id should not be empty",
         call. = FALSE)
  }
  
  # read the content of files paths provided in ths column
  output <- list()
  output$samples <- list()
  for (i in 1:nrow(input_table
                    )) {
    
    # show progress bar indicating how many samples out of total were processed
    if (verbose) {
      message(paste0("Processing sample ",i," out of ",nrow(input_table)))
    }
    
    output$samples[[input_table$sample_id[i]]] <- list(data = read_polya_single(input_table$polya_path[i], ...))
    output$samples[[input_table$sample_id[i]]]$meta <- input_table[i, -which(names(input_table) %in% c("polya_path"))]
  }
  
  
  output$metadata_table <- input_table
  
    if (process_references) {
    # get mapped references for all samples and store in the "references" element of the list
    if (verbose) {
      message("Processing references")
    }
    output$references <- data.frame(reference=get_references(output,reference_column = "reference"),symbol=NA)
    
    # convert all reference_columns in data to factors
    for (i in 1:length(output$samples)) {
      output$samples[[i]]$data$reference <- as.factor(output$samples[[i]]$data$reference)
    }
  }
  
  if (verbose) {
    message("Finished all samples from the input table")
  }
  
  # set additional class for the output, inicating nanotail object (list) with poly(A) data
  class(output) <- c(class(output), "nanotail_polya_data")
  
  return(output)
}


