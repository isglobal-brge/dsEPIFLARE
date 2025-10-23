#' Check if required methods exist on server
#'
#' @param client dsHPC client object
#' @param required_methods Character vector of required method names
#' @return List with status and missing methods
#' @keywords internal
check_server_methods <- function(client, required_methods) {
  tryCatch({
    # Get available methods from server
    methods <- client$getMethods()
    
    # getMethods returns a list where names are method names
    available_methods <- names(methods)
    
    # Check which required methods are missing
    missing <- setdiff(required_methods, available_methods)
    
    if (length(missing) > 0) {
      return(list(
        status = FALSE,
        available = available_methods,
        missing = missing,
        message = paste0("Missing methods on server: ", paste(missing, collapse = ", "))
      ))
    }
    
    return(list(
      status = TRUE,
      available = available_methods,
      missing = character(0),
      message = "All required methods available"
    ))
    
  }, error = function(e) {
    return(list(
      status = FALSE,
      available = character(0),
      missing = required_methods,
      message = paste("Error checking methods:", e$message)
    ))
  })
}

#' Validate resource object
#'
#' @param resource Resource object with connection information
#' @return TRUE if valid, error otherwise
#' @keywords internal
validate_resource <- function(resource) {
  required_fields <- c("name", "url", "format", "secret")
  
  if (!is.list(resource)) {
    stop("resource must be a list")
  }
  
  missing_fields <- setdiff(required_fields, names(resource))
  if (length(missing_fields) > 0) {
    stop("resource missing required fields: ", paste(missing_fields, collapse = ", "))
  }
  
  if (resource$format != "dshpc.api") {
    stop("resource format must be 'dshpc.api', got: ", resource$format)
  }
  
  return(TRUE)
}

#' NULL-coalescing operator
#'
#' @param x First value
#' @param y Default value if x is NULL
#' @return x if not NULL, otherwise y
#' @keywords internal
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

