
#' Plot dose response curve
#'
#' Plot dose response curve fit using fit_drc_model function.
#' @param model 4-parameter drc model.
#' @param dr_data Dose response dataframe. If columns are not labelled "DOSE" and
#'  "RESPONSE", indicate the indices of the relevant columns.
#' @param dose_col Index of dose/concentration column.
#' @param resp_col Index of response column. % Viability recommended.
#' @param ic50 IC50 value for highlight on plot.
#'
#' @return Plot of dose response relationship.
#' @export
#'
#' @examples
plot_dr_curve <- function(model,
                          dr_data,
                          dose_col = NULL,
                          resp_col = NULL,
                          ic50 = NULL){

  if (!all("DOSE" %in% colnames(dr_data), "RESPONSE" %in% colnames(dr_data))) {
    if (is.null(dose_col) | is.null(resp_col)){
      return("Please enter dose and response column indices.")
    } else {
      colnames(dr_data)[dose_col] <- "DOSE"
      colnames(dr_data)[resp_col] <- "RESPONSE"
    }
  }

  newdata<-expand.grid(conc=exp(seq(log(min(dr_data$DOSE)),
                                    log(max(dr_data$DOSE)),
                                    length=100))) # data for smooth curve
  curve_bulked <- predict(model, newdata= newdata, interval = "confidence")
  newdata$p <- curve_bulked[,1]
  newdata$pmax <- curve_bulked[,2]
  newdata$pmin <- curve_bulked[,3]

  plot <- ggplot(dr_data,aes(x=DOSE,y=RESPONSE))+
    geom_point()+
    geom_ribbon(data=newdata,
                aes(x=conc,y=p,ymin=pmin,ymax=pmax), stat="identity",alpha=0.2)+
    geom_line(data=newdata,
              aes(x=conc,y=p))+
    scale_x_log10()+
    ylab("Viability(% Control)")+
    ylim(c(0,120))+
    geom_hline(yintercept=50, linetype="dashed",
               color = "red", size=1)+
    {if(!is.null(ic50))geom_vline(xintercept = ic50, linetype="dashed",
                                  color = "blue", size=1)} +
    theme_bw()

  return(plot)
}


