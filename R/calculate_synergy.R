calculate_cdi <- function(dr_dataset,
                          drug1_col,
                          drug2_col){

  dr_dataset$cdi <- dr_dataset$co_treatment / dr_dataset[drug1_col] * dr_dataset[drug2_col]
  dr_dataset$log2cdi <- log(dr_dataset$cdi)
  return(dr_dataset)
}


calculate_bliss <- function(dr_dataset,
                            drug1_col,
                            drug2_col){

  bliss <- dr_dataset[drug1_col] + dr_dataset[drug2_col] - dr_dataset[drug1_col] * dr_dataset[drug2_col]
  dr_dataset$bliss_independence <- bliss
  return(dr_dataset)
}
