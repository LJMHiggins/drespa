#' Fit dose response curve
#'
#' @description
#' Function to fit 4 parameter dose response model.
#'
#' @param dr_data Dose response dataframe. If columns are not labelled "DOSE" and
#'  "RESPONSE", indicate the indices of the relevant columns.
#' @param dose_col Index of dose/concentration column.
#' @param resp_col Index of response column. Percentage viability recommended.
#'
#' @return A 4-parameter drc model if convergence is achieved. Otherwise
#'  returns NULL.
#' @export
#'
#' @examples
fit_drc_model <- function(dr_data,
                          dose_col = 1,
                          resp_col = 2){
  # House keeping
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

compute_abs_IC50 <- function(){

}

#' @title Calculate AUC or AAC
#'
#' @description
#' Calculate the Area Under the Curve (AUC) or the Area Above the Curve (AAC) for a dose-response curve.
#'
#' @param dr_data Dose-response dataframe. If columns are not labeled "DOSE" and "RESPONSE," indicate the indices of the relevant columns.
#' @param dose_col Index of the dose/concentration column.
#' @param resp_col Index of the response column, percentage viability recommended.
#' @param aac Logical. Whether or not to return the Area Above the Dose-Response Curve (AAC).
#'
#' @return
#' A numeric value representing the AUC or AAC.
#'
#' @examples
#' \dontrun{
#'   data <- data.frame(
#'     DOSE = c(0, 1, 2, 3, 4),
#'     RESPONSE = c(100, 90, 80, 70, 60)
#'   )
#'   computeAUC(data, dose_col = 1, resp_col = 2)
#' }
#'
#' @export
calculateAUC <- function(dr_data,
                         dose_col,
                         resp_col,
                         aac = FALSE) {
  x <- dr_data[[dose_col]]
  y <- dr_data[[resp_col]]
  if (length(x) != length(y)) {
    print("Impute vectors must be of same length.")
  }
  idx <- order(x)
  x <- x[idx]
  y <- y[idx]
  conc_diff <- diff(x)
  areas <-  sum((rowMeans(cbind(y[-length(y)], y[-1]))) * conc_diff)
  max_auc <- sum(100 * conc_diff)
  normalized_auc <- areas / max_auc

  if (aac){
    return(1 - normalized_auc)
  } else {
    return(normalized_auc)
  }
}

#' Calculate Drug Sensitivity Scores including DSS1, DSS2 and DSS3
#'
#' @param model
#' @param df
#'
#' @return
#' @export
#'
#' @examples
computeDSS <- function(model, df){
  ## Needs adapting
  data = data.frame(
    ic50 = model$coefficients[[4]],
    slope = model$coefficients[[1]],
    MIN_resp = min(df$RESPONSE),
    MAX_resp = max(df$RESPONSE),
    Cmin = min(df$DOSE),
    Cmax = max(df$DOSE)
  )
  dss_1 <- DSS::DSS(data, concn.scale = 1e-6, log = TRUE, DSS.type = 1)
  dss_2 <- DSS::DSS(data, concn.scale = 1e-6, log = TRUE, DSS.type = 2)
  dss_3 <- DSS::DSS(data, concn.scale = 1e-6, log = TRUE, DSS.type = 3)

  results <- data.frame(DSS_1 = dss_1,
                        DSS_2 = dss_2,
                        DSS_3 = dss_3)
  colnames(results) <- c("DSS_1", "DSS_2", "DSS_3")

  return(results)
}


fit_glm_model <- function(dr_data,
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
  # House keeping
}




