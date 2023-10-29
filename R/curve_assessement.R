#' Curve assessment - QC metrics and DSS calculation
#'
#' Extract curve parameter estimates along with QC metrics including R2, RMSE
#'  and parameter confidence intervals extracted from drc model parameters. This
#'  includes EC_50 (half of maximum inhibition - also known as relative IC50).
#'  Additionally, Drug Sensitivity Scores are calculated (DSS1, DSS2 & DSS3). For this,
#'  ensure data is percentage of viability (100 percent = DMSO control) on a scale of 0 - 100,
#'  and that concentration is in uM and not in log10 scale as transformation is
#'  applied within function.
#'
#' @param model 4-parameter drc model.
#' @param dr_data Dose response dataframe. If columns are not labelled "DOSE" and
#'  "RESPONSE", indicate the indices of the relevant columns.
#' @param dose_col Index of dose/concentration column.
#' @param resp_col Index of response column. Percentage viability recommended.
#'
#' @return Dataframe with single row.
#' @export
#'
#' @examples
run_curve_assessment <- function(model,
                                 dr_data,
                                 dose_col = NULL,
                                 resp_col = NULL){

  if (!all("DOSE" %in% colnames(dr_data), "RESPONSE" %in% colnames(dr_data))) {
    if (is.null(dose_col) | is.null(resp_col)){
      return("Please enter dose and response column indices.")
    } else {
      colnames(dr_data)[dose_col] <- "DOSE"
      colnames(dr_data)[resp_col] <- "RESPONSE"
    }
  }
  #### R-squared value
  .calculate.r.sq <- function(model, df){
    observed <- df$RESPONSE
    predicted <- fitted(model)
    residuals <- observed - predicted
    SSR <- sum(residuals^2)
    SST <- sum((observed - mean(observed))^2)
    R_squared <- 1 - (SSR / SST)

    return(R_squared)
  }
  #### RMSE function
  .calculate.rmse <- function(model, df){
    observed <- df$RESPONSE
    predicted <- fitted(model)
    residuals <- observed - predicted
    RMSE <- sqrt(mean(residuals^2))

    return(RMSE)
  }

  #### Confidence interval for parameters
  .parameter.estimates <- function(model){
    # For a 4 param curve fit
    df <- broom::tidy(model, conf.int = TRUE) #Extracts relative IC50 and other params in table format
    df$CI_range <- abs(df$conf.low - df$conf.high)
    ## Return as single row ##
    df.final <- df[-2]
    df.final %>%
      reshape2::melt()%>%
      unite(param_id, term, variable, sep = "_")%>%
      pivot_wider(names_from = param_id, values_from = value) -> out

    return(out)
  }
  #### Calculate ABSOLUTE IC50 (also known as GIC50)
  .calculate_gIC50 <- function(model){
    concentrations <- exp(seq(log(min(model$data$DOSE)), log(max(model$data$DOSE)), length.out = 1000))
    predicted_responses <- predict(model, newdata = data.frame(DOSE = concentrations))
    # Get concentration where the predicted response is closest to 50%
    index_closest_to_50 <- which.min(abs(predicted_responses - 50))
    gIC50 <- concentrations[index_closest_to_50]

    return(gIC50)
  }
  #### Calculate DSS scores (1 - 3)
  .calculate.DSS <- function(model, df){
    data = data.frame(
      ic50 = model$coefficients[[4]], # Relative IC50 (or EC50)
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

  #### Combine and return results from all sub functions
  r_sq <- .calculate.r.sq(model = model, df = dr_data)
  rmse <- .calculate.rmse(model = model, df = dr_data)
  gIC50 <- .calculate_gIC50(model = model)
  metrics <- .parameter.estimates(model = model)
  metrics$r_squared <- r_sq
  metrics$rmse <- rmse
  metrics$gIC50 <- gIC50

  dss_scores <- .calculate.DSS(model = model, df = dr_data)
  # could insert calculation AAC here also
  metrics <- cbind(metrics, dss_scores)
  return(metrics)
}
