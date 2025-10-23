#' Read and Process IDAT Files
#'
#' @description
#' Processes methylation array IDAT files using the dsHPC infrastructure.
#' The function uploads a tar.gz archive containing IDAT files and phenotype data,
#' then processes them on the remote server using minfi.
#'
#' @param resource A list containing connection information with fields:
#'   \itemize{
#'     \item name: Resource name (character)
#'     \item url: Server URL (character)
#'     \item format: Must be "dshpc.api" (character)
#'     \item secret: API key for authentication (character)
#'   }
#' @param idat_object Raw vector containing the tar.gz file data (use readBin to load file)
#' @param filename Original filename of the tar.gz file (for server reference)
#' @param output_format Output format: "hybrid" (Parquet+JSON, default), "json", or "parquet"
#'
#' @return A list containing:
#'   \itemize{
#'     \item status: "success" or "error"
#'     \item data: List with betas (DataFrame) and pheno (DataFrame)
#'     \item summary: Processing statistics
#'     \item format: Data format used
#'   }
#'
#' @details
#' The input tar.gz file must contain:
#' \itemize{
#'   \item pheno.csv: Phenotype data with columns PID, Sex, array_id, etc.
#'   \item IDATs/: Folder with IDAT files
#' }
#'
#' The function automatically:
#' \itemize{
#'   \item Checks if read_idat method exists on server
#'   \item Uploads the file (with deduplication)
#'   \item Processes IDAT files with minfi
#'   \item Applies functional normalization
#'   \item Filters samples and CpGs based on detection p-values
#'   \item Returns beta values and phenotype data
#' }
#'
#' @examples
#' \dontrun{
#' # Define resource
#' resource <- list(
#'   name = "epiflare_hpc",
#'   url = "http://localhost:8001",
#'   format = "dshpc.api",
#'   secret = "your_api_key_here"
#' )
#'
#' # Read file into memory
#' idat_file <- "path/to/your_data.tar.gz"
#' idat_data <- readBin(idat_file, "raw", file.info(idat_file)$size)
#' 
#' # Process IDAT files
#' result <- ds.readIDAT(
#'   resource = resource,
#'   idat_object = idat_data,
#'   filename = basename(idat_file),
#'   output_format = "hybrid"
#' )
#'
#' # Access results
#' betas <- result$data$betas
#' pheno <- result$data$pheno
#' print(result$summary)
#' }
#'
#' @export
#' @import dsHPC
ds.readIDAT <- function(resource, 
                        idat_object, 
                        filename,
                        output_format = "hybrid") {
  
  # Validate inputs
  validate_resource(resource)
  
  if (!is.raw(idat_object)) {
    stop("idat_object must be a raw vector (use readBin to load file)")
  }
  
  if (length(idat_object) == 0) {
    stop("idat_object is empty")
  }
  
  if (missing(filename) || is.null(filename) || filename == "") {
    stop("filename must be provided")
  }
  
  if (!output_format %in% c("hybrid", "json", "parquet")) {
    stop("output_format must be one of: hybrid, json, parquet")
  }
  
  # Initialize client
  message("Initializing dsHPC client...")
  client <- HPCResourceClient$new(resource)
  
  # Check if read_idat method exists
  message("Checking server methods...")
  method_check <- check_server_methods(client, "read_idat")
  
  if (!method_check$status) {
    stop(method_check$message)
  }
  
  message("✓ Method 'read_idat' available on server")
  
  # Create API config
  # Extract host and port from URL
  url_parts <- strsplit(resource$url, ":")[[1]]
  base_url <- paste(url_parts[1:2], collapse = ":")  # http://host
  port <- as.integer(url_parts[3])
  
  api_config <- dsHPC:::create_api_config(
    base_url = base_url,
    port = port,
    api_key = resource$secret,
    auth_header = "X-API-Key",
    auth_prefix = ""
  )
  
  # Get file info
  file_size_mb <- length(idat_object) / 1024 / 1024
  message(sprintf("Processing %s (%.2f MB)...", filename, file_size_mb))
  
  # Submit job
  message("Submitting job to server...")
  
  method_chain <- list(
    list(
      method_name = "read_idat",
      parameters = list(
        output_format = output_format
      )
    )
  )
  
  start_time <- Sys.time()
  
  result <- execute_processing_chain(
    config = api_config,
    content = idat_object,
    method_chain = method_chain,
    upload_filename = filename
  )
  
  end_time <- Sys.time()
  duration <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  # Check for errors
  if (!is.null(result$status) && result$status == "error") {
    stop("Job failed: ", result$error %||% "Unknown error")
  }
  
  message(sprintf("✓ Processing completed in %.2f seconds", duration))
  
  # Add timing information
  result$timing <- list(
    duration_seconds = duration,
    duration_minutes = duration / 60
  )
  
  return(result)
}

