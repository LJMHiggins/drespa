#' Plot random dose response curves
#'
#' Plot dose response curves from large dataset using random selection. Intended
#'  for checking large datasets with many dose response relationships too large
#'  to plot all possible curves. For example, if user wants to perform visual
#'  checks on model fits, this function could be used to generate a sample of
#'  dose response graphs construction.
#'  Note to self, in the future include option for different model types.
#'
#' @param cell_line_col Index for cell line column (if not named "CELL_LINE_NAME").
#' @param drug_name_col Index for drug name column (if not named "DRUG_NAME").
#' @param dose_col Index of dose/concentration column.
#' @param resp_col Index of response column.
#' @param sample_size Number of dose response curves to plot. Note this does not
#'  account for lack of convergence.
#' @param save_pdf Logical value, determines whether to save output pdf.
#' @param report_title String to be used as pdf file name.
#' @param dr_dataset Data frame with dose response data in long format (Column for
#'  cell line id, drug id, dose and percentage viability).
#'
#' @return
#' @export
#'
#' @examples
plot_random_drs <- function(cell_line_col = NULL,
                            drug_name_col = NULL,
                            dose_col = NULL,
                            resp_col = NULL,
                            sample_size = 25,
                            save_pdf = TRUE,
                            report_title = "Random_curves",
                            dr_dataset){
  ## Generic col checker ##
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
      # Avoid issue with duplicate columns
      dplyr::select(
        dr_dataset,
        c(CELL_LINE_NAME, DRUG_NAME, DOSE, RESPONSE)) -> dr_dataset
    }
  }
  ## Random df subset ##
  grid_exec <- expand.grid(unique(dr_dataset$CELL_LINE_NAME),
                           unique(dr_dataset$DRUG_NAME),
                           stringsAsFactors = FALSE)

  grid_random <- grid_exec[sample(nrow(grid_exec), sample_size), ]

  plot_list <- list()

  for (i in 1:nrow(grid_random)){
    df_sub <- subset(dr_dataset, CELL_LINE_NAME == grid_random[i,1] & DRUG_NAME == grid_random[i,2])
    skip_to_next <- FALSE

    tryCatch({
      model <- drespa::fit_drc_model(dr_data = df_sub,
                                     dose_col = dose_col,
                                     resp_col = resp_col)
      plot <- drespa::plot_dr_curve(model = model,
                                    dr_data = df_sub,
                                    dose_col = dose_col,
                                    resp_col = resp_col)
      plot <- plot + ggtitle(paste0("C, D: ", grid_random[i,1], ", ", grid_random[i,2]))

      plot_list[[paste0(grid_random[i, 1], grid_random[i, 2])]] <- plot
    }, error = function(e) {skip_to_next <<- TRUE}
    )
    if(skip_to_next){next} # In case of model fit fail
  }
  if (save_pdf){
    jumps <- seq(from =1, to = length(plot_list), by=9)
    pdf(paste0(report_title, ".pdf"))
    for (i in 1:length(jumps)){
      # Check if less than 9 plots are left
      if (jumps[i]+8 > length(plot_list)){
        plots_to_save <- plot_list[jumps[i]: length(plot_list)]
        plot(cowplot::plot_grid(plotlist = plots_to_save, label_size = 11))
        break
      }
      plots_to_save <- plot_list[jumps[i]: (jumps[i]+8)]
      plot(cowplot::plot_grid(plotlist = plots_to_save, label_size = 11))
    }
    dev.off()
  }

  return(plot_list)
}
