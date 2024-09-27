#' Convert a list with poly(A) predictions to long format data.frame with metadata columns (like in the original nanoTail)
#'
#' @description
#' When data are loaded by the Nanotail package thery are stored as list of class "nanotail_polya_data", 
#' which is convenient for internal processing but not for further analysis, including visualization. 
#' This function converts the nanotail_polya_data list to a data.frame, containing all per-read poly(A) data
#' with metadata stored as separate columns.
#' 
#' @param input_list a list - output of read_polya_multiple() with poly(A) predictions
#'
#' @return data.frame (tibble)
#' @export
#' 
#' @examples
#' \dontrun{
#' 
#' polya_list_to_data_frame(input_list)
#' 
#' }
#'
polya_list_to_data_frame <- function(input_list) {
  # check if input_list is provided
  if (missing(input_list)) {
    stop("List is missing. Please provide a valid list argument",
         call. = FALSE)
  }
  
  # check if input_list is a nanotail class
  if (!is.nanotail_polya_data(input_list)) {
    stop("Input should be provided as a nanotail_polya_data list. Please load your data with read_polya_multiple()",
         call. = FALSE)
  }
  
  
  
  # remove elements from samples, if data is empty (nrow == 0)
  input_list$samples <- lapply(input_list$samples, function(x) {
    if (nrow(x$data) == 0) {
      return(NULL)
    } else {
      return(x)
    }
  })
  
  # add metadata to data
  for (i in 1:length(input_list$samples)) {
    input_list$samples[[i]]$data <- cbind(input_list$samples[[i]]$meta, input_list$samples[[i]]$data)
  }
  
  # return a data.frame with all elements of the list containg metadata columns
  output <- do.call(rbind, lapply(input_list$samples, function(x) x$data))
  row.names(output) <- NULL # get rid of row names
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
get_transcript_df_from_polya_list <- function(input_list, transcript_ids,transcript_id_column="transcript",...) {
  # check if input_list is provided
  if (missing(input_list)) {
    stop("List is missing. Please provide a valid list argument",
         call. = FALSE)
  }
  # check if input_list is a nanotail class
  if (!is.nanotail_polya_data(input_list)) {
    stop("Input should be provided as a nanotail_polya_data list. Please load your data with read_polya_multiple()",
         call. = FALSE)
  }
  
  # check if transcript is provided
  if (missing(transcript_ids)) {
    stop("Transcript is missing. Please provide a valid transcript argument",
         call. = FALSE)
  }
  
  
  # check if transcript is a character vector or single character
  if (!is.character(transcript_ids) & length(transcript) != 1) {
    stop("Transcript should be a character vector or single character",
         call. = FALSE)
  }
  
  

  # filter a data element of each element of the list, to keep only rows where transcript column is equal to the transcript argument
  # transcript column is specified by the transcript_id_column argument

  input_list <- filter_polya_list_by_transcript(input_list = input_list,transcript_ids = transcript_ids,transcript_id_column = transcript_id_column,...)
  
  output_df <- polya_list_to_data_frame(input_list,...)
  
  return(output_df)
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

summarize_polya_list <- function(input_list,transcript_ids=NA,transcript_id_column="transcript",...) {
  # check if input_list is provided
  if (missing(input_list)) {
    stop("List is missing. Please provide a valid list argument",
         call. = FALSE)
  }
  
  # check if input_list is a nanotail class
  if (!is.nanotail_polya_data(input_list)) {
    stop("Input should be provided as a nanotail_polya_data list. Please load your data with read_polya_multiple()",
         call. = FALSE)
  }

 
  
  if (!is.na(transcript_ids)) {
    message(paste0("Filtering transcripts by ",transcript_id_column," = ",transcript_ids))
    input_list <- filter_polya_list_by_transcript(input_list,transcript_id_column=transcript_id_column,transcript_ids=transcript_ids,...)
  }
  
  # calculate - number of elements in each data element, mean, median, sd, min, max of polya_length column
  output <- lapply(input_list$samples, function(x) {
    
    data <- x$data
    meta <- x$meta
    n <- nrow(data)
    mean_polya <- mean(data$polya_length)
    median_polya <- median(data$polya_length)
    sd_polya <- sd(data$polya_length)
    min_polya <- min(data$polya_length)
    max_polya <- max(data$polya_length)
    if (!is.na(transcript_ids)) { # add transcript names to output, if was specified in function call.
      transcript_ids <- paste0(transcript_ids,sep=",")
      meta[[transcript_id_column]] <- transcript_ids
    }
    output <- c(n, mean_polya, median_polya, sd_polya, min_polya, max_polya)
    names(output) <- c("n", "mean_polya", "median_polya", "sd_polya", "min_polya", "max_polya")
    output <- c(meta, output)
    return(output)
  })
  
  # add content of metadata for each sample as first columns of the output
  output <- do.call(rbind, output)
  rownames(output) <- NULL
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
  # check if input_list is a nanotail class
  if (!is.nanotail_polya_data(input_list)) {
    stop("Input should be provided as a nanotail_polya_data list. Please load your data with read_polya_multiple()",
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
  
 
  
  # filter metadata columns of each element of the list
  output_samples <- lapply(input_list$samples, function(x) {
    meta <- x$meta
    meta <- meta[metadata]
    return(list(data = x$data, meta = meta))
  })
  
  input_list$samples <- output_samples
  input_list$metadata_table <- input_list$metadata_table[metadata]
  
  return(input_list)
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
filter_polya_list_by_transcript <- function(input_list, transcript_ids,transcript_id_column="transcript",verbose=TRUE) {
  # check if input_list is provided
  if (missing(input_list)) {
    stop("List is missing. Please provide a valid list argument",
         call. = FALSE)
  }
  # check if input_list is a nanotail class
  if (!is.nanotail_polya_data(input_list)) {
    stop("Input should be provided as a nanotail_polya_data list. Please load your data with read_polya_multiple()",
         call. = FALSE)
  }
  
  # check if transcript is provided
  if (missing(transcript_ids)) {
    stop("Transcript is missing. Please provide a valid transcript argument",
         call. = FALSE)
  }
  
  # check if transcript is a character vector or single character
  if (!is.character(transcript_ids) & length(transcript_ids) != 1) {
    stop("Transcript should be a character vector or single character",
         call. = FALSE)
  }
  
  # get references from references data.frame based on provided transcript ids
  if (transcript_id_column=='reference') {
    filtered_references <- transcript_ids
    if (verbose) {
      message("Using provided reference ids for data filtering")
    }
  }
  else {
    references_table <- input_list$references 
    filtered_references <- references_table[references_table[[transcript_id_column]] %in% transcript_ids,]$reference
    if (verbose) {
      message(paste0("Got reference ids using provided transcript_ids stored in column",transcript_id_column))
    }
  }
  
  output_samples <- lapply(input_list$samples, function(x) {
    data <- x$data
    meta <- x$meta
    data <- data[data[["reference"]] %in% filtered_references,]
    return(list(data = data, meta = meta))
    if (verbose) {
      message(paste0("Filtered data for ",length(transcript_ids)," transcript ids"))
    }
  })
  
  input_list$samples <- output_samples
  return(input_list)

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
  # check if input_list is a nanotail class
  if (!is.nanotail_polya_data(input_list)) {
    stop("Input should be provided as a nanotail_polya_data list. Please load your data with read_polya_multiple()",
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
  
  output <- unique(unlist(lapply(input_list$samples, function(x) x$data[[reference_column]])))
  return(output)
}



#' Filter samples from the list of poly(A) predictions
#' 
#' @param input_list a list - output of read_polya_multiple() with poly(A) predictions
#' @param samples a character vector with sample names to keep
#' 
#' @return list with filtered samples
#' @export
#' 
#' @examples
#' \dontrun{
#' 
#' filter_polya_list_samples(input_list,samples=c("sample1","sample2"))
#' 
#' }
#' 
filter_polya_list_samples <- function(input_list, samples) {
  
  # check if inout_list is of nanotail_polya_data class
  # if (!is.nanotail_polya_data(input_list)) {
  #   stop("Input should be provided as a nanotail_polya_data class",
  #        call. = FALSE)
  # }
  
  output_list <- input_list
  
  #if samples are provided, leave only the listed samples in tje output list
  if (!missing(samples)) {
    message("samples")
    if (!is.character(samples)) {
      stop("Samples should be provided as a character vector",
           call. = FALSE)
    }
    
    # check if each element of samples is present in meta column sample_name
    if (!all(samples %in% input_list$samples[[1]]$meta$sample_name)) {
      stop("Some of the samples are not present in the input list",
           call. = FALSE)
    }
    
    message("Filtering")
    #filter the list, leaving only samples which are listed in the meta
    output_samples <- lapply(input_list$samples, function(x) {
      meta <- x$meta
      
      if (meta$sample_name %in% samples) {
        message(meta$sample_name)
        return(x)
      }
      else {
        message(paste("filter out",meta$sample_name))
        return(NULL)
      }
    })
    
    #remove all NULL elements from output_samples list
    output_samples <- output_samples[sapply(output_samples, function(x) !is.null(x))]
    
    output_list$samples <- output_samples
  
  }  
  return(output_list)
}

# filter_list is a list of pairs of column names and values from the metadata element of input_list
# values should be provided as vectors
get_samples_list_from_metadata <- function(input_list,filter_list) {
  
  # check if inout_list is of nanotail_polya_data class
  # if (!is.nanotail_polya_data(input_list)) {
  #   stop("Input should be provided as a nanotail_polya_data class",
  #        call. = FALSE)
  # }
  
  # check if filter_list is provided
  if (missing(filter_list)) {
    stop("Filter list is missing. Please provide a valid filter_list argument",
         call. = FALSE)
  }
  
  # check if filter_list is a list
  if (!is.list(filter_list)) {
    stop("Filter list should be provided as a list",
         call. = FALSE)
  }
  
  # # check if each element of the list has another list named meta
  # if (!all(sapply(filter_list, function(x) length(x) == 2))) {
  #   stop("Each element of the list should have a pair of column name and value",
  #        call. = FALSE)
  # }
  
  # # check if each element of the list has another list named meta
  # if (!all(sapply(filter_list, function(x) is.character(x[[1]])))) {
  #   stop("First element of the pair should be a character vector",
  #        call. = FALSE)
  # }
  # 
  # # check if each element of the list has another list named meta
  # if (!all(sapply(filter_list, function(x) is.character(x[[2]])))) {
  #   stop("Second element of the pair should be a character vector",
  #        call. = FALSE)
  # }
  # 
  # check if each element of the list has another list named meta
  if (!names(filter_list) %in% names(input_list$samples[[1]]$meta)) {
    stop("First element of the pair should be a column name from the metadata",
         call. = FALSE)
  }
  
  
  # Using provided filters, return a character vector of samples names (column sample_name in meta) that match all of the criteria,
  # we are not interested in the whole samples element, with data, just the vector of name
  output <- lapply(input_list$samples, function(x) {
    meta <- x$meta
    for (i in 1:length(filter_list)) {
      meta <- meta[meta[[names(filter_list)[[i]]]] %in% filter_list[[i]],]
    }
    return(meta$sample_name)
  })
  
  

  return(output)
}

# functions which merges two nanotail_polya_data lists
merge_nanotail_polya_data <- function(input1,input2,verbose=FALSE) {
  
  # check if input1 is provided
  if (missing(input1)) {
    stop("First list is missing. Please provide a valid list argument",
         call. = FALSE)
  }
  
  # check if input1 is a nanotail class
  if (!is.nanotail_polya_data(input1)) {
    stop("First input should be provided as a nanotail_polya_data list. Please load your data with read_polya_multiple()",
         call. = FALSE)
  }
  
  if(!validate_nanotail_polya_data(input1)) {
    stop("First input is not a valid nanotail_polya_data list. Please load your data with read_polya_multiple()",
         call. = FALSE)
  }
  
  # check if input2 is provided
  if (missing(input2)) {
    stop("Second list is missing. Please provide a valid list argument",
         call. = FALSE)
  }
  
  # check if input2 is a nanotail class
  if (!is.nanotail_polya_data(input2)) {
    stop("Second input should be provided as a nanotail_polya_data list. Please load your data with read_polya_multiple()",
         call. = FALSE)
  }
  
  if(!validate_nanotail_polya_data(input2)) {
    stop("Second input is not a valid nanotail_polya_data list. Please load your data with read_polya_multiple()",
         call. = FALSE)
  }

  
  # check if input1 and input2 metadata_tables has the same columns in the same order
  # if not, add missing columns to metadata table, and store information about what is missing for each input
  # if the same, just return variable identical_metadata=TRUE
  identical_metadata=TRUE
  
  add_missing_columns <- function(input,missing_columns) {
    for (i in 1:length(missing_columns)) {
      input[[missing_columns[[i]]]] <- NA
    }
    return(input)
  }
  
  if (!identical(names(input1$metadata_table),names(input2$metadata_table))) {
    if (verbose) {
      message("Metadata tables have different columns. Creating matching metadata columns")
    }
    identical_metadata=FALSE
    missing_columns_input2 <- setdiff(names(input1$metadata_table),names(input2$metadata_table))
    if (length(missing_columns_input2) > 0) {
      
      for (i in 1:length(missing_columns_input2)) {
        input2$metadata_table[[missing_columns_input2[[i]]]] <- NA
      }
      
      for (i in 1:length(input2$samples)) {
        input2$samples[[i]]$meta <- add_missing_columns(input2$samples[[i]]$meta,missing_columns_input2)
      }
      
    }
    missing_columns_input1 <- setdiff(names(input2$metadata_table),names(input1$metadata_table))
    if (length(missing_columns_input1) > 0) {
      
      for (i in 1:length(missing_columns_input1)) {
        input1$metadata_table[[missing_columns_input1[[i]]]] <- NA
      }
      # add missing column also to each meta element in the input1
      for (i in 1:length(input1$samples)) {
        input1$samples[[i]]$meta <- add_missing_columns(input1$samples[[i]]$meta,missing_columns_input1)
      }
    }
  }
  
  # Reorganize both metadata_tables from input1 and input2, so they have same order of columns
  
  if (!identical_metadata) {
    input2$metadata_table <- input2$metadata_table[names(input1$metadata_table)]
  }

  # check if input1 and input2 have the same references objects
  # if the same - just take the references from input1
  # if not the same, merge two reference data.frames, leaving only unique rows
  
  if (identical(input1$references,input2$references)) {
    output_references <- input1$references
  }
  else {
    if (verbose) {
      message("References are different. Merging")
    }
    output_references <- rbind(input1$references,input2$references)
    output_references <- unique(output_references)
  }
  
  
  # check if there are element (by name) in input2 which are also present in input1
  # If names are duplicated, check also the data element from both inputs for identity
  # if the same, skip second sample, if different, assign different sample_id (suffix dupl) and issue the warning
  # remember to use the new sample_id in the meta element
  
  if (verbose) {
    message("Merging samples...")
  }
  if (any(names(input2$samples) %in% names(input1$samples))) {
    for (i in 1:length(names(input2$samples))) {
      if (names(input2$samples)[[i]] %in% names(input1$samples)) {
        if (identical(input1$samples[[names(input2$samples)[[i]]]]$data,input2$samples[[names(input2$samples)[[i]]]]$data)) {
          if (verbose) {
            message(paste0("Sample ",names(input2$samples)[[i]]," is present in both inputs. Leaving only the first one."))
          }
          # remove duplicated sample from input2 metadata_table, by sample_id
          input2$metadata_table <- input2$metadata_table[!input2$metadata_table$sample_id %in% names(input2$samples)[[i]],]
          next
        }
        else {
         
          new_name <- paste0(names(input2$samples)[[i]],"_dupl")
          #input2$samples[[new_name]] <- input2$samples[[names(input2$samples)[[i]]]]
          if (verbose) {
            message(paste0("Sample ",names(input2$samples)[[i]]," is already present in the first list, but the data are different. Assigning new sample_id: ",new_name))
          }
          input2$metadata_table[input2$metadata_table$sample_id %in% names(input2$samples)[[i]],]$sample_id <- new_name
          names(input2$samples)[[i]] <- new_name
          input2$samples[[i]]$meta$sample_id <- new_name
          # assign also the new name duplicated sample from input2 metadata_table, by sample_id
          
        }
      }
    }
  }
  
  
  
  # merge samples
  output_samples <- c(input1$samples,input2$samples)

  
  # merge metadata_tables, removing duplicated samples
  output_metadata_table <- rbind(input1$metadata_table,input2$metadata_table)
  
  
  # return the merged list
  
  output <- input1
  output$samples <- output_samples
  output$metadata_table <- output_metadata_table
  output$references <- output_references
  
  if (verbose) {
    message("Merging completed")
  }
  
  return(output)
  
  
}
