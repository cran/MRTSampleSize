#' @title  Calculate power for micro-randomized trials
#'
#' @description This function calculates power for micro-randomized trials (MRTs)
#' based on methodology developed in Sample Size Calculations for Micro-randomized Trials in
#' mHealth by Liao et al. (2016) <DOI:10.1002/sim.6847>.
#'
#' @param days Duration of the study.
#' @param occ_per_day Number of decision time points per day.
#' @param prob Randomization probability, i.e. the probability of assigning the treatment at a decision time
#'             point. This can be constant, or time-varying probabilities can be specified by a vector
#'             specifying randomization probabilities for each day or decision time.
#' @param beta_shape The trend for the proximal treatment effect, choices are constant, linear or quadratic.
#' Note:
#'  \enumerate{
#'              \item{Constant} The proximal treatment effect stays constant over the study.
#'              \item{Linear} The linearly increasing form of a proximal treatment effect might be
#'              used if participants will get more enthusiastically engage in the apps and thus
#'              the proximal effect will increase as the study goes. The linearly decreasing form
#'              of a proximal treatment effect might be used if participants are likely to disengage
#'              the activity suggestionss and thus the proximal effect will decrease as the study goes.
#'              \item{Quadratic}  The quadratic form of a proximal treatment effect might be used if
#'              you expect that initially participants will enthusiastically engage in the apps and
#'              thus the proximal effect will get higher. Then, as the study goes on, some participants
#'              are likely to disengage or begin to ignore the activity suggestions and hence a downward
#'              trend.
#'  }
#'
#' @param beta_mean Average of proximal treatment effect.
#' @param beta_initial Initial value of proximal treatment effect when beta_shape is linear or quadratic.
#' @param beta_quadratic_max Day of maximal proximal treatment effect when beta_shape is quadratic.
#' @param tau_shape The pattern for expected availability, choices are constant, linear or quadratic.
#' Note:
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
#'
#' @param tau_mean Average of expected availability.
#' @param tau_initial Initial Value of expected availability when tau_shape is linear or quadratic.
#' @param tau_quadratic_max Changing point of availability when tau_shape is quadratic.
#' @param sample_size Number of participants
#' @param sigLev Significance level
#'
#' @return power.
#'
#' @export calculatePower
#' @references Seewald, N.J.; Sun, J.; Liao, P. "MRT-SS Calculator: An R Shiny Application for Sample Size
#' Calculation in Micro-Randomized Trials". arXiv:1609.00695
#' @examples
#'     calculatePower(days=42,
#'                    occ_per_day=5,
#'                    prob=0.4,
#'                    beta_shape="quadratic",
#'                    beta_mean=0.1,
#'                    beta_initial=0,
#'                    beta_quadratic_max=28,
#'                    tau_shape="quadratic",
#'                    tau_mean=0.5,
#'                    tau_initial=0.7,
#'                    tau_quadratic_max=42,
#'                    sample_size=40,
#'                    sigLev=0.05)
#'
#'     prob1 <- c(replicate(35,0.7),replicate(35,0.6),replicate(35,0.5),replicate(35,0.4))
#'     calculatePower(days=28,
#'                    occ_per_day=5,
#'                    prob=prob1,
#'                    beta_shape="quadratic",
#'                    beta_mean=0.1,
#'                    beta_initial=0,
#'                    beta_quadratic_max=28,
#'                    tau_shape="quadratic",
#'                    tau_mean=0.5,
#'                    tau_initial=0.7,
#'                    tau_quadratic_max=42,
#'                    sample_size=40,
#'                    sigLev=0.05)#'
#'

calculatePower <- function(days,occ_per_day,prob,
                           beta_shape,beta_mean,beta_initial,beta_quadratic_max,
                           tau_shape,tau_mean,tau_initial,tau_quadratic_max,
                           sample_size,sigLev){

  if(sample_size != round(sample_size)){
    stop("Error: Please enter integer value for Number of Participants")
  }
  if(sample_size <= 0){
    stop("Error: Please specify Number of Participants greater than 0")
  }

  validateParameters(days,occ_per_day,prob,
                     beta_shape,beta_mean,beta_initial,beta_quadratic_max,
                     tau_shape,tau_mean,tau_initial,tau_quadratic_max,
                     sigLev)


  beta_input <- generateBeta(days,occ_per_day,beta_shape,beta_mean,beta_initial,beta_quadratic_max)
  tau_input <- generateTau(days,occ_per_day,tau_shape,tau_mean,tau_initial,tau_quadratic_max)

  ### calculate power ###

  total <- days*occ_per_day
  input_avail <- vector('numeric', total)
  input_effect <- vector('numeric',total)

  p_input <- 3
  if(beta_shape == "constant"){
    p_input <- 1
  } else if(beta_shape == "linear"){
    p_input <- 2
  }
  ### We assume that the proximal treatment effect is consant on each day ###
  for(k in 1:days)
  {
    input_effect[(occ_per_day*k-occ_per_day+1):(occ_per_day*k)] = replicate(occ_per_day, beta_input[k])
  }

  ### We assume that the expected availability is constant on each day ###
  for(k in 1:days)
  {
    input_avail[(occ_per_day*k-occ_per_day+1):(occ_per_day*k)] = replicate(occ_per_day, tau_input[k])
  }

  if(length(prob) == 1 || length(prob) == days || length(prob) == days*occ_per_day){

    power <- PowerCalculation(days, occ_per_day, sample_size, input_effect, input_avail, delta = prob, alpha0=sigLev, p=p_input, q=3)

    if(power >= 0.5){
      print(sprintf("The power is %5.2f%s with sample size %3d when the significance level is %4.2f.", power*100,"%",sample_size,sigLev))
    }else{
      ### if the calculated power is less than 50%, output an warning ###
      print(sprintf("The power is less than 50%s with sample size %3d when the significance level is %4.2f.","%", sample_size, sigLev))
    }

  }else{
    print(paste("Parameter Prob has wrong length!"))
  }
  return(power)

}
