# ========================== PIPELINE LOG SINK =================================
#' @name log_sink
#' @description saves the messages and outputs from the console in a text file
#' useful to keep logs/track during parallel processing
#' @param FILE a file connection in form of file("name", open = "a|wt") where
#' "a" stands for "append" and "wt" writes a new file.
#' @param START if TRUE, opens the sink, if FALSE, closes the sink

log_sink <- function(FILE, START){
  if(START == TRUE){
    # --- 1. Start the sink
    sinkfile <- FILE
    sink(sinkfile, type = "output", append = TRUE)
    sink(sinkfile, type = "message", append = TRUE)
    return(sinkfile)
  } else {
    # --- 2. Close connection and stop sink
    sink(type = "message")
    sink(type = "output")
    close(FILE)
  }
} # END FUNCTION
