#' Convert a list with poly(A) predictions to long format data.frame with metadata columns (like in the original nanoTail)
#'
#' @param input_list a list - output of read_polya_multiple() with poly(A) predictions
#'
#' @return data.frame (tibble)
#' @export
#'
polya_list_to_data_frame <- function(input_list) {
  # check if input_list is provided
  if (missing(input_list)) {
    stop("List is missing. Please provide a valid list argument",
         call. = FALSE)
  }
  # check if input_list is a list
  if (!is.list(input_list)) {
    stop("List should be provided as a list",
         call. = FALSE)
  }
  

  # check if each element of the list has another list named meta
  if (!all(sapply(input_list, function(x) "meta" %in% names(x)))) {
    stop("Each element of the list should have another list named meta",
         call. = FALSE)
  }
  
  # check if each element of the list has another data element
  if (!all(sapply(input_list, function(x) "data" %in% names(x)))) {
    stop("Each element of the list should have another data element",
         call. = FALSE)
  }
  
  # check if each element of the list has a data element which is a data.frame
  if (!all(sapply(input_list, function(x) is.data.frame(x$data)))) {
    stop("Each element of the list should have a data element which is a data.frame",
         call. = FALSE)
  }
  
  # check if each element of the list has a meta element which is a list
  if (!all(sapply(input_list, function(x) is.list(x$meta)))) {
    stop("Each element of the list should have a meta element which is a list",
         call. = FALSE)
  }
  
  # skip empty data.frames
  input_list <- input_list[sapply(input_list, function(x) nrow(x$data) > 0)]
  
  # add metadata to data
  for (i in 1:length(input_list)) {
    input_list[[i]]$data <- cbind(input_list[[i]]$meta, input_list[[i]]$data)
  }
  
  # return a data.frame with all elements of the list containg metadata columns
  output <- do.call(rbind, lapply(input_list, function(x) x$data), make.row.names = FALSE)
  return(output)
}


#' Get data for a specific transcripts from a list of poly(A) predictions
#'
#' @param input_list a list - output of read_polya_multiple() with poly(A) predictions
#' @param transcript a character vector or single character with transcript names to filter
#'
#' @return data.frame (tibble) with data for selected transcripts
#' @export
#'
#' @examples
get_transcript_data_from_polya_list <- function(input_list, transcript,transcript_id_column="transcript") {
  # check if input_list is provided
  if (missing(input_list)) {
    stop("List is missing. Please provide a valid list argument",
         call. = FALSE)
  }
  # check if input_list is a list
  if (!is.list(input_list)) {
    stop("List should be provided as a list",
         call. = FALSE)
  }
  
  # check if transcript is provided
  if (missing(transcript)) {
    stop("Transcript is missing. Please provide a valid transcript argument",
         call. = FALSE)
  }
  
  # check if transcript is a character vector or single character
  if (!is.character(transcript) & length(transcript) != 1) {
    stop("Transcript should be a character vector or single character",
         call. = FALSE)
  }
  
  
  # check if each element of the list has another list named meta
  if (!all(sapply(input_list, function(x) "meta" %in% names(x)))) {
    stop("Each element of the list should have another list named meta",
         call. = FALSE)
  }
  
  # check if each element of the list has another data element
  if (!all(sapply(input_list, function(x) "data" %in% names(x)))) {
    stop("Each element of the list should have another data element",
         call. = FALSE)
  }
  
  # check if each element of the list has a data element which is a data.frame
  if (!all(sapply(input_list, function(x) is.data.frame(x$data)))) {
    stop("Each element of the list should have a data element which is a data.frame",
         call. = FALSE)
  }
  
  # check if each element of the list has a meta element which is a list
  if (!all(sapply(input_list, function(x) is.list(x$meta)))) {
    stop("Each element of the list should have a meta element which is a list",
         call. = FALSE)
  }
  
  # filter a data element of each element of the list, to keep only rows where transcript column is equal to the transcript argument
  # transcript column is specified by the transcript_id_column argument

  output <- lapply(input_list, function(x) {
    data <- x$data
    meta <- x$meta
    data <- data[data[[transcript_id_column]] %in% transcript,]
    return(list(data = data, meta = meta))
  })
  
  
  
  # add metadata to data
  
  for (i in 1:length(output)) {
    output[[i]]$data <- cbind(output[[i]]$meta, output[[i]]$data,)
  }
  
  # return a data.frame
  message("output data.frame")
  output <- do.call(rbind, lapply(input_list, function(x) x$data))
  return(output)
  
}


#' Summarize poly(A) predictions from a list
#' 
#' @param input_list a list - output of read_polya_multiple() with poly(A) predictions
#'  
#' 
#' @return data.frame (tibble) with summary statistics for each sample
#' @export
#' 
#' @examples
#' \dontrun{
#' 
#' summarize_polya_list(input_list)
#' 
#' }
#' 

summarize_polya_list <- function(input_list,transcript=NA,transcript_id_column="transcript") {
  # check if input_list is provided
  if (missing(input_list)) {
    stop("List is missing. Please provide a valid list argument",
         call. = FALSE)
  }
  # check if input_list is a list
  if (!is.list(input_list)) {
    stop("List should be provided as a list",
         call. = FALSE)
  }

  # check if each element of the list has another list named meta
  if (!all(sapply(input_list, function(x) "meta" %in% names(x)))) {
    stop("Each element of the list should have another list named meta",
         call. = FALSE)
  }
  
  # check if each element of the list has another data element
  if (!all(sapply(input_list, function(x) "data" %in% names(x)))) {
    stop("Each element of the list should have another data element",
         call. = FALSE)
  }
  
  # check if each element of the list has a data element which is a data.frame
  if (!all(sapply(input_list, function(x) is.data.frame(x$data)))) {
    stop("Each element of the list should have a data element which is a data.frame",
         call. = FALSE)
  }
  
  # check if each element of the list has a meta element which is a list
  if (!all(sapply(input_list, function(x) is.list(x$meta)))) {
    stop("Each element of the list should have a meta element which is a list",
         call. = FALSE)
  }
  
  # calculate - number of elements in each data element, mean, median, sd, min, max of polya_length column
  output <- lapply(input_list, function(x) {
    
    if (!is.na(transcript)) {
      message("Filtering transcript",transcript)
      data <- x$data
      data <- data[data[[transcript_id_column]] %in% transcript,]
      x$data <- data
    }
    else {
      data <- x$data
    }
    meta <- x$meta
    n <- nrow(data)
    mean_polya <- mean(data$polya_length)
    median_polya <- median(data$polya_length)
    sd_polya <- sd(data$polya_length)
    min_polya <- min(data$polya_length)
    max_polya <- max(data$polya_length)
    if (!is.na(transcript)) { # add transcript name to output, if was specified in function call.
      meta$transcript <- transcript
    }
    output <- c(n, mean_polya, median_polya, sd_polya, min_polya, max_polya)
    names(output) <- c("n", "mean_polya", "median_polya", "sd_polya", "min_polya", "max_polya")
    output <- c(meta, output)
    return(output)
  })
  
  # add content of metadata for each sample as first columns of the output
  output <- do.call(rbind, output)
  return(output)
}

#' Filter metadata from the list of poly(A) predictions
#' 
#' @param input_list a list - output of read_polya_multiple() with poly(A) predictions
#' @param metadata a character vector with names of metadata columns to keep
#' 
#' @return list with filtered metadata
#' @export
#' 
#' @examples
#' \dontrun{
#' 
#' filter_metadata(input_list,metadata=c("sample_name","group"))
#' 
#' }
#' 
drop_polya_list_metadata <- function(input_list, metadata) {
  # check if input_list is provided
  if (missing(input_list)) {
    stop("List is missing. Please provide a valid list argument",
         call. = FALSE)
  }
  # check if input_list is a list
  if (!is.list(input_list)) {
    stop("List should be provided as a list",
         call. = FALSE)
  }
  
  # check if metadata is provided
  if (missing(metadata)) {
    stop("Metadata is missing. Please provide a valid metadata argument",
         call. = FALSE)
  }
  
  # check if metadata is a character vector
  if (!is.character(metadata)) {
    stop("Metadata should be a character vector",
         call. = FALSE)
  }
  
  # check if each element of the list has another list named meta
  if (!all(sapply(input_list, function(x) "meta" %in% names(x)))) {
    stop("Each element of the list should have another list named meta",
         call. = FALSE)
  }
  
  # check if each element of the list has another data element
  if (!all(sapply(input_list, function(x) "data" %in% names(x)))) {
    stop("Each element of the list should have another data element",
         call. = FALSE)
  }
  
  # check if each element of the list has a data element which is a data.frame
  if (!all(sapply(input_list, function(x) is.data.frame(x$data)))) {
    stop("Each element of the list should have a data element which is a data.frame",
         call. = FALSE)
  }
  
  # check if each element of the list has a meta element which is a list
  if (!all(sapply(input_list, function(x) is.list(x$meta)))) {
    stop("Each element of the list should have a meta element which is a list",
         call. = FALSE)
  }
  
  # filter metadata columns of each element of the list
  output <- lapply(input_list, function(x) {
    meta <- x$meta
    meta <- meta[metadata]
    return(list(data = x$data, meta = meta))
  })
  
  return(output)
}




#' Filter data from the list of poly(A) predictions
#' 
#' @param input_list a list - output of read_polya_multiple() with poly(A) predictions
#' @param transcript a character vector or single character with transcript names to filter
#' @param transcript_id_column a character vector with column name with transcript names
#' 
#' @return list with filtered data
#' @export
#' 
#' @examples
#' \dontrun{
#' 
#' filter_data(input_list,transcript=c("ACTB"),transcript_id_column="transcript")
#' 
#' }
filter_polya_list_by_transcript <- function(input_list, transcript,transcript_id_column="transcript") {
  # check if input_list is provided
  if (missing(input_list)) {
    stop("List is missing. Please provide a valid list argument",
         call. = FALSE)
  }
  # check if input_list is a list
  if (!is.list(input_list)) {
    stop("List should be provided as a list",
         call. = FALSE)
  }
  
  # check if transcript is provided
  if (missing(transcript)) {
    stop("Transcript is missing. Please provide a valid transcript argument",
         call. = FALSE)
  }
  
  # check if transcript is a character vector or single character
  if (!is.character(transcript) & length(transcript) != 1) {
    stop("Transcript should be a character vector or single character",
         call. = FALSE)
  }
  
  
  # check if each element of the list has another list named meta
  if (!all(sapply(input_list, function(x) "meta" %in% names(x)))) {
    stop("Each element of the list should have another list named meta",
         call. = FALSE)
  }
  
  # check if each element of the list has another data element
  if (!all(sapply(input_list, function(x) "data" %in% names(x)))) {
    stop("Each element of the list should have another data element",
         call. = FALSE)
  }
  
  # check if each element of the list has a data element which is a data.frame
  if (!all(sapply(input_list, function(x) is.data.frame(x$data)))) {
    stop("Each element of the list should have a data element which is a data.frame",
         call. = FALSE)
  }
  
  # check if each element of the list has a meta element which is a list
  if (!all(sapply(input_list, function(x) is.list(x$meta)))) {
    stop("Each element of the list should have a meta element which is a list",
         call. = FALSE)
  }
  
  # filter a data element of each element of the list, to keep only rows where transcript column is equal to the transcript argument
  # transcript column is specified by the transcript_id_column argument
  
  output <- lapply(input_list, function(x) {
    data <- x$data
    meta <- x$meta
    data <- data[data[[transcript_id_column]] %in% transcript,]
    return(list(data = data, meta = meta))
  })
  
  return(output)
  
}

#' Get references from the list of poly(A) predictions
#' 
#' @param input_list a list - output of read_polya_multiple() with poly(A) predictions
#' @param reference_column a character vector with column name with references
#' 
#' @return list with references
#' @export
#' 
#' @examples
#' \dontrun{
#' 
#' get_references(input_list,reference_column="reference")
#' 
#' }
#'  
get_references <- function(input_list,reference_column="reference") {

  # check if input_list is provided
  if (missing(input_list)) {
    stop("List is missing. Please provide a valid list argument",
         call. = FALSE)
  }
  # check if input_list is a list
  if (!is.list(input_list)) {
    stop("List should be provided as a list",
         call. = FALSE)
  }

  # check if each element of the list has another list named meta
  if (!all(sapply(input_list, function(x) "meta" %in% names(x)))) {
    stop("Each element of the list should have another list named meta",
         call. = FALSE)
  }
  
  # check if each element of the list has another data element
  if (!all(sapply(input_list, function(x) "data" %in% names(x)))) {
    stop("Each element of the list should have another data element",
         call. = FALSE)
  }
  
  # check if each element of the list has a data element which is a data.frame
  if (!all(sapply(input_list, function(x) is.data.frame(x$data)))) {
    stop("Each element of the list should have a data element which is a data.frame",
         call. = FALSE)
  }
  
  # check if each element of the list has a meta element which is a list
  if (!all(sapply(input_list, function(x) is.list(x$meta)))) {
    stop("Each element of the list should have a meta element which is a list",
         call. = FALSE)
  }
  
  # check if transcript is provided
  if (missing(reference_column)) {
    stop("Reference_column is missing. Please provide a valid reference_column argument",
         call. = FALSE)
  }
  
  # check if transcript is a string, not the vector
  
  if (!is.character(reference_column) & length(reference_column) != 1) {
    stop("Reference_column should be a character vector or single character",
         call. = FALSE)
  }
  
  
  # get the content of reference column of data element of each element of the list
  # combine all the references for each list elements into single vector, with all values unique
  # return the produced vector
  
  output <- unique(unlist(lapply(input_list, function(x) x$data[[reference_column]])))
  
}

