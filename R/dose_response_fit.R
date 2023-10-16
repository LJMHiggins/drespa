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


#' Calculate AUC
#'
#' @param dr_data Dose response dataframe. If columns are not labelled "DOSE" and
#'  "RESPONSE", indicate the indices of the relevant columns.
#' @param dose_col Index of dose/concentration column.
#' @param resp_col Index of response column. % Viability recommended.
#' @param aac Boolean. Whether or not to return area above the dr curve.
#'
#' @return AUC or AAC.
#' @export
#'
#' @examples
computeAUC <- function(dr_data,
                    dose_col = 1,
                    resp_col = 2,
                    aac = FALSE) {

  curve_data <- data.frame(viab=dr_data[1], conc=dr_data[1],)
  sorted_curve <- curve_data[order(curve_data$conc),]

  viab_diff <- diff(sorted_data$viab)
  conc_diff <- diff(sorted_data$conc)

  areas <- 0.5 * (viab_diff[1:(length(viab_diff) - 1)] + viab_diff[2:length(viab_diff)]) * conc_diff
  area_total <- sum(areas)
  normalized_auc <- area_total / log(sorted_data[nrow(sorted_data), ]$conc / sorted_data[1, ]$conc)

  if (aac){
    return(1 - normalized_auc)
  } else {
    return(normalized_auc)
  }
}

computeDSS <- function(model, df){
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




