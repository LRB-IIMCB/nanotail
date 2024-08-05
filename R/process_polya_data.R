# having a list, which has elements with names of samples, inside each element there is another list named meta, where each element of this list describes the metadata of the sample
# another element for each sample name is called data, and is a data.frame wit numeric data, in 8 columns
# write a function, which for each sample in a list takes the metadata, and adds them as additional columns to data data.frame
# return a data.frame with all elements of the list containg metadata columns
# use base R function, dplyr or purrr usage is prohibited

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
  
  # add metadata to data
  for (i in 1:length(input_list)) {
    input_list[[i]]$data <- cbind(input_list[[i]]$meta, input_list[[i]]$data,)
  }
  
  # return a data.frame with all elements of the list containg metadata columns
  output <- do.call(rbind, lapply(input_list, function(x) x$data))
  return(output)
}


# write a function which will take two arguments - input_list and transcript
# a function will filter a data element of each element of the list, to keep only rows where transcript column is equal to the elements of vector provided as a transcript argument
# additional columns should be added to the output, taken from the metadata list - element of list for each sample
# return a data.frame
# use base R function, dplyr or purrr usage is prohibited

#' Get data for a specific transcripts from a list of poly(A) predictions
#'
#' @param input_list a list - output of read_polya_multiple() with poly(A) predictions
#' @param transcript a character vector or single character with transcript names to filter
#'
#' @return data.frame (tibble) with data for selected transcripts
#' @export
#'
#' @examples
get_transcript_data_from_list <- function(input_list, transcript) {
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
  for (i in 1:length(input_list)) {
    input_list[[i]]$data <- input_list[[i]]$data[input_list[[i]]$data$transcript %in% transcript, ]
  }
  
  # add metadata to data
  for (i in 1:length
       (input_list)) {
    input_list[[i]]$data <- cbind(input_list[[i]]$meta, input_list[[i]]$data)
  }
  
  # return a data.frame
  output <- do.call(rbind, lapply(input_list, function(x) x$data))
  return(output)
  
}


# return summary of each sample in the input list
# calculate - numbr of elements in each data element, mean, median, sd, min, max of polya_length column
# add content of metadata for each sample as first columns of the output
# use base R function, dplyr or purrr usage is prohibited


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

summarize_polya_list <- function(input_list) {
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
  
  # calculate - numbr of elements in each data element, mean, median, sd, min, max of polya_length column
  output <- lapply(input_list, function(x) {
    data <- x$data
    meta <- x$meta
    n <- nrow(data)
    mean_polya <- mean(data$polya_length)
    median_polya <- median(data$polya_length)
    sd_polya <- sd(data$polya_length)
    min_polya <- min(data$polya_length)
    max_polya <- max(data$polya_length)
    output <- c(n, mean_polya, median_polya, sd_polya, min_polya, max_polya)
    names(output) <- c("n", "mean_polya", "median_polya", "sd_polya", "min_polya", "max_polya")
    output <- c(meta, output)
    return(output)
  })
  
  # add content of metadata for each sample as first columns of the output
  output <- do.call(rbind, output)
  return(output)
}





