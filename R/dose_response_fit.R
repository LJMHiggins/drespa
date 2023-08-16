
#' Fit dose response curve
#'
#' Function to fit 4 parameter dose response model.
#'
#' @param dr_data Dose response dataframe. If columns are not labelled "DOSE" and
#'  "RESPONSE", indicate the indices of the relevant columns.
#' @param dose_col Index of dose/concentration column.
#' @param resp_col Index of response column. % Viability recommended.
#'
#' @return A 4-parameter drc model if convergence is achieved. Otherwise
#'  returns NULL.
#' @export
#'
#' @examples
fit_drc_model <- function(dr_data,
                          dose_col = 1,
                          resp_col = 2){

  if (!all("DOSE" %in% colnames(dr_data), "RESPONSE" %in% colnames(dr_data))) {
    if (is.null(dose_col) | is.null(resp_col)){
      return("Please enter dose and response column indices.")
    } else {
      colnames(dr_data)[dose_col] <- "DOSE"
      colnames(dr_data)[resp_col] <- "RESPONSE"
    }
  }

  output <- tryCatch(
    {
      model <- drc::drm(formula = RESPONSE ~ DOSE,
                        data = dr_data,
                        fct=drc::LL.4(names = c("hill", "min_value", "max_value", "IC_50")))
    },
    error = function(e){
      print("Could not fit curve.")
      return(NULL)
    }
  )
  return(output)
}
