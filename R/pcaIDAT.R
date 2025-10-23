#' Read IDAT Files and Perform PCA Analysis
#'
#' @description
#' Processes methylation array IDAT files and performs Principal Component Analysis (PCA)
#' using the dsHPC infrastructure. The heavy computation (read_idat + PCA calculation) 
#' happens on the server, while plot generation occurs locally for flexibility.
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
#' @param top_CpG Number of most variable CpGs to use for PCA (default: 10000)
#' @param array Array type: "450k", "EPIC", or "EPICv2" (default: "EPIC")
#' @param plot_title Custom title for the PCA plot (default: "PCA Plot")
#' @param color_by Metadata variable to color points by (default: "Case_Cont")
#' @param shape_by Metadata variable to shape points by (default: "Sex")
#'
#' @return A list containing:
#'   \itemize{
#'     \item plot: ggplot2 plot object
#'     \item pca_coordinates: DataFrame with PCA coordinates for all samples
#'     \item sample_metadata: DataFrame with sample metadata
#'     \item variance_explained: Named vector with variance explained by each PC
#'     \item summary: Processing statistics (n_samples, n_cpgs_used, etc.)
#'     \item timing: Execution time information
#'   }
#'
#' @details
#' The function executes a processing chain on the server:
#' \enumerate{
#'   \item read_idat: Processes IDAT files with minfi
#'   \item pca_analysis: Performs PCA on beta values
#' }
#'
#' The server returns:
#' \itemize{
#'   \item PCA coordinates (all PCs)
#'   \item Sample metadata
#'   \item Variance explained by each PC
#' }
#'
#' Plot generation happens locally, allowing for customization without reprocessing data.
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
#' # Process IDAT files and perform PCA
#' result <- ds.pcaIDAT(
#'   resource = resource,
#'   idat_object = idat_data,
#'   filename = basename(idat_file),
#'   top_CpG = 10000,
#'   array = "EPIC",
#'   plot_title = "My PCA Analysis"
#' )
#'
#' # Display plot
#' print(result$plot)
#'
#' # Save plot
#' ggsave("pca_plot.png", result$plot, width = 10, height = 8)
#'
#' # Access PCA coordinates
#' pca_coords <- result$pca_coordinates
#'
#' # Check variance explained
#' print(result$variance_explained[1:5])
#'
#' # Create custom plot with different PCs
#' library(ggplot2)
#' plot_data <- cbind(result$sample_metadata, result$pca_coordinates)
#' ggplot(plot_data, aes(x = PC3, y = PC4, color = Condition)) +
#'   geom_point()
#' }
#'
#' @export
#' @import dsHPC
#' @import ggplot2
#' @importFrom rlang .data
ds.pcaIDAT <- function(resource,
                       idat_object,
                       filename,
                       top_CpG = 10000,
                       array = "EPIC",
                       plot_title = "PCA Plot",
                       color_by = "Case_Cont",
                       shape_by = "Sex") {
  
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
  
  if (!array %in% c("450k", "EPIC", "EPICv2")) {
    stop("array must be one of: 450k, EPIC, EPICv2")
  }
  
  # Initialize client
  message("Initializing dsHPC client...")
  client <- HPCResourceClient$new(resource)
  
  # Check if both methods exist
  message("Checking server methods...")
  required_methods <- c("read_idat", "pca_analysis")
  method_check <- check_server_methods(client, required_methods)
  
  if (!method_check$status) {
    stop(method_check$message)
  }
  
  message(sprintf("✓ Methods available: %s", paste(required_methods, collapse = ", ")))
  
  # Create API config
  url_parts <- strsplit(resource$url, ":")[[1]]
  base_url <- paste(url_parts[1:2], collapse = ":")
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
  
  # Define processing chain
  message("Submitting processing chain to server...")
  message("  Step 1: read_idat (process IDAT files)")
  message("  Step 2: pca_analysis (compute PCA)")
  
  method_chain <- list(
    list(
      method_name = "read_idat",
      parameters = list(
        output_format = "hybrid"  # Use Parquet for efficiency
      )
    ),
    list(
      method_name = "pca_analysis",
      parameters = list(
        top_CpG = as.integer(top_CpG),
        array = array
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
  
  server_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  
  # Check for errors first
  if (!is.null(result$status) && result$status == "error") {
    stop("Job failed: ", result$error %||% "Unknown error")
  }
  
  # Check if we have data
  if (is.null(result$data)) {
    stop("No data returned from server")
  }
  
  message(sprintf("✓ Server processing completed in %.2f seconds", server_time))
  
  # Extract data from result
  message("Generating plot locally...")
  
  # Get variance explained
  var_exp <- result$data$variance_explained
  if (is.list(var_exp)) {
    var_exp <- unlist(var_exp)
  }
  
  # Combine PCA coordinates with metadata for plotting
  pca_coords <- as.data.frame(result$data$pca_coordinates)
  metadata <- as.data.frame(result$data$sample_metadata)
  plot_data <- cbind(metadata, pca_coords)
  
  # Check if requested variables exist
  if (!color_by %in% names(metadata)) {
    warning(sprintf("Variable '%s' not found in metadata. Available: %s", 
                    color_by, paste(names(metadata), collapse = ", ")))
    color_by <- NULL
  }
  
  if (!shape_by %in% names(metadata)) {
    warning(sprintf("Variable '%s' not found in metadata. Available: %s",
                    shape_by, paste(names(metadata), collapse = ", ")))
    shape_by <- NULL
  }
  
  # Create plot
  plot_start <- Sys.time()
  
  plt <- ggplot(plot_data, aes(x = PC1, y = PC2))
  
  if (!is.null(color_by) && !is.null(shape_by)) {
    plt <- plt + geom_point(aes(col = .data[[color_by]], shape = .data[[shape_by]]), 
                           size = 2.5, alpha = 0.9)
    plt <- plt + stat_ellipse(aes(col = .data[[color_by]]), 
                              level = 0.99, alpha = 0.5, lty = 2, show.legend = FALSE)
  } else if (!is.null(color_by)) {
    plt <- plt + geom_point(aes(col = .data[[color_by]]), size = 2.5, alpha = 0.9)
    plt <- plt + stat_ellipse(aes(col = .data[[color_by]]),
                              level = 0.99, alpha = 0.5, lty = 2, show.legend = FALSE)
  } else if (!is.null(shape_by)) {
    plt <- plt + geom_point(aes(shape = .data[[shape_by]]), size = 2.5, alpha = 0.9)
  } else {
    plt <- plt + geom_point(size = 2.5, alpha = 0.9)
  }
  
  plt <- plt + 
    labs(
      title = plot_title,
      subtitle = sprintf("%d samples, %s CpGs", 
                        result$data$n_samples,
                        format(result$data$n_cpgs_used, big.mark = ",")),
      x = sprintf("PC1 (%.2f%% variance)", var_exp["PC1"]),
      y = sprintf("PC2 (%.2f%% variance)", var_exp["PC2"])
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "gray40")
    )
  
  plot_time <- as.numeric(difftime(Sys.time(), plot_start, units = "secs"))
  
  message(sprintf("✓ Plot generated in %.3f seconds", plot_time))
  
  # Return comprehensive result
  return(list(
    plot = plt,
    pca_coordinates = pca_coords,
    sample_metadata = metadata,
    variance_explained = var_exp,
    summary = list(
      n_samples = result$data$n_samples,
      n_cpgs_used = result$data$n_cpgs_used,
      n_cpgs_filtered = result$data$n_cpgs_filtered
    ),
    timing = list(
      server_seconds = server_time,
      server_minutes = server_time / 60,
      plot_seconds = plot_time,
      total_seconds = server_time + plot_time
    )
  ))
}

