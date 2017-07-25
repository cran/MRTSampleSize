#' @title  plot the graph for the expected availability
#'
#' @description plot of the graphs for the expected availability, i.e., the expected
#' probability that a participant is available to receive treatment at a decision time.
#' when the pattern for the expected availability is constant, linear or quadractic.
#'
#' @param days Duration of the study.
#' @param occ_per_day Number of decision time points per day.
#' @param tau_shape The pattern for expected availability, choices are constant, linear or quadratic.
#'  Note:
#'  \enumerate{
#'              \item{Constant} The expected availability stays constant over the study.
#'              \item{Linear} A linearly increasing pattern of expected availability might be used if
#'              participants will find the intervention useful and thus more likely to turn the
#'              intervention on. A linearly decreasing pattern of expected availability might be used
#'              if participants learn more about the intervetion and get bored through the course of
#'              the study and thus getting less likely to turn on the invervention.
#'              \item{Quadratic} A quadratic pattern of availability. Here the changing point of
#'              availability refers to day of either maximal of minimal availability, depending on the
#'              input values of initial and average availability.
#'  }
#' @param tau_mean Average of expected availability.
#' @param tau_initial Initial Value of expected availability when tau_shape is linear or quadratic.
#' @param tau_quadratic_max Changing point of availability when tau_shape is quadratic.
#'
#' @importFrom graphics plot abline legend
#' @return A graph for expected availability.
#'
#' @export plotExpectAvail
#' @examples
#'    plotExpectAvail(days=42,
#'                    occ_per_day=5,
#'                    tau_shape="quadratic",
#'                    tau_mean=0.5,
#'                    tau_initial=0.7,
#'                    tau_quadratic_max=42)
#'
#'
plotExpectAvail <- function(days,occ_per_day,tau_shape,tau_mean,tau_initial,tau_quadratic_max){
  if(days != round(days)){
    stop("Error: Please enter integer values for the number of days")
  }
  if(days <= 0){
    stop("Error: Please specify the number of days greater than 0")
  }
  if(occ_per_day != round(occ_per_day)){
    stop("Error:Please enter integer for the number of occasions per day")
  }
  if(occ_per_day <= 0){
    stop("Error: Please specify the number of occasions per day greater than 0")
  }
  validateEffectAvailParameters(tau_shape,tau_mean)

  tau_shape < tolower(tau_shape)
  tau_input <- generateTau(days,occ_per_day,tau_shape,tau_mean,tau_initial,tau_quadratic_max)
  if(tau_shape == "constant"){
    plot(tau_input,xlab = "Days", ylab = "Availability", ylim = c(0,1), type = "o",
         pch = 16, cex = 0.8, col = 4)
    abline(h=tau_mean, lty = 2)
    legend("topleft", legend=c('Availability','Average Availability'), col = c(4,1),lty = c(1,2), pch=c(16,NA),bty = "n")
  } else if(tau_shape == "linear"){

    if(min(tau_input) < 0){
      warning("Warning: Some values of the availability are less than 0.")
    }
    if(max(tau_input) > 1){
      warning("Warning: Some values of the availability are bigger than 1.")
    }
    plot(tau_input,xlab = "Days", ylab = "Availability", ylim = c(0,1), type = "o",
         pch = 16, cex = 0.8, col = 4)
    abline(h=tau_mean, lty = 2)
    legend("topleft", legend=c('Availability','Average of Availability'), col = c(4,1),lty = c(1,2), pch=c(16,NA),bty = "n")

  } else if(tau_shape == "quadratic"){
    if(min(tau_input) < 0){
      warning("Warning: Some values of the availability are less than 0.")
    }
    if(max(tau_input) > 1){
      warning("Warning: Some values of the availability are bigger than 1.")
    }
    plot(tau_input, xlab = "Days", ylab = "Availability", ylim = c(0,1), type = "o",
         pch = 16, cex = 0.8, col = 4)
    abline(h=tau_mean, lty = 2)
    legend("topleft", legend=c('Availability','Average of Availability'), col = c(4,1),lty = c(1,2), pch=c(16,NA),bty = "n")
  }

}
