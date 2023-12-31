#' Parallel processing application of curve assessment function
#'
#' Wrapper function for large-scale implementation of model fitting and
#' curve assessment.
#'
#' @param cell_line_col Index for the cell line column (if not named "CELL_LINE_NAME").
#' @param drug_name_col Index for the drug name column (if not named "DRUG_NAME").
#' @param dr_dataset Data frame with dose-response data in long format (Columns for
#'   cell line id, drug id, dose, and percentage viability).
#' @param dose_col Index of the dose/concentration column.
#' @param resp_col Index of the response column.
#' @param n_cores Number of CPUs to use.
#' @param output_fname File name for CSV file output. If NULL, the file will not be
#'   saved by the function.
#'
#' @return Data frame with curve metrics and DSS scores. Each row pertains to a
#'   single dose curve (identified by cell line/drug id columns in the output).
#'
#' @examples
#' \dontrun{
#' curve_process_pipeline(
#'   cell_line_col = 1,
#'   drug_name_col = 2,
#'   dose_col = 3,
#'   resp_col = 4,
#'   dr_dataset = your_data_frame,
#'   n_cores = 4,
#'   output_fname = "output.csv"
#' )
#' }
#'
#' @export
curve_process_pipeline <- function(cell_line_col = NULL,
                                   drug_name_col = NULL,
                                   dose_col = NULL,
                                   resp_col = NULL,
                                   dr_dataset,
                                   n_cores = 4,
                                   output_fname = NULL) {

  library(parallel)
  library(doParallel)

  if (!all("CELL_LINE_NAME" %in% colnames(dr_dataset),
           "DRUG_NAME" %in% colnames(dr_dataset),
           "DOSE" %in% colnames(dr_dataset),
           "RESPONSE" %in% colnames(dr_dataset))) {
    if (is.null(cell_line_col) | is.null(drug_name_col) | is.null(dose_col) | is.null(resp_col)){
      return("Please enter cell line name, drug name, dose, and response column indices.")
    } else {
      colnames(dr_dataset)[cell_line_col] <- "CELL_LINE_NAME"
      colnames(dr_dataset)[drug_name_col] <- "DRUG_NAME"
      colnames(dr_dataset)[dose_col] <- "DOSE"
      colnames(dr_dataset)[resp_col] <- "RESPONSE"
    }
  }

  grid_exec <- expand.grid(unique(dr_dataset$CELL_LINE_NAME),
                           unique(dr_dataset$DRUG_NAME),
                           stringsAsFactors = FALSE)

  n_cores <- n_cores
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  clusterEvalQ(cl, {
    library(drc)
    library(dplyr)
    library(tidyr)
  })

  start_time <- Sys.time()

  df_results <- foreach::foreach(i= 1:nrow(grid_exec), .combine = "rbind") %dopar% {

    df_sub <- subset(dr_dataset, CELL_LINE_NAME == grid_exec[i,1] & DRUG_NAME == grid_exec[i,2])

    model <- drespa::fit_drc_model(dr_data = df_sub,
                                   dose_col = dose_col,
                                   resp_col = resp_col)

    if (!is.null(model)){
      metrics <- drespa::run_curve_assessment(model = model,
                                              dr_data = df_sub,
                                              dose_col = dose_col,
                                              resp_col = resp_col)
      metrics$CELL_LINE_NAME <- grid_exec[i,1]
      metrics$DRUG_NAME <- grid_exec[i,2]
      metrics
    }
  }
  if (!is.null(output_fname)) {
    write.csv(df_results, output_fname)
  }

  end_time <- Sys.time()
  stopCluster(cl)
  registerDoSEQ()

  time_taken <- end_time - start_time
  print(paste0("Complete! Process took: ", time_taken))

  return(df_results)
}
