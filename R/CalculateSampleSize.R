#' @title  Calculate sample size for micro-randomized trials
#'
#' @description This function calculates the sample size for micro-randomized trials (MRTs)
#' based on methodology developed in Sample Size Calculations for Micro-randomized Trials in
#' mHealth by Liao et al. (2016) <DOI:10.1002/sim.6847>.
#'
#' @param days The duration of the study.
#' @param occ_per_day The number of decision time points per day.
#' @param prob The randomization probability, i.e. the probability of assigning the treatment at a decision time
#'             point. This can be constant, or time-varying probabilities can be specified by by a vector
#'             specifying randomization probabilities for each day or decision time.
#' @param beta_shape The trend for the proximal treatment effect; choices are constant, linear or quadratic.
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
#' @param beta_mean The average of proximal treatment effect.
#' @param beta_initial The initial value of proximal treatment effect when beta_shape is linear or quadratic.
#' @param beta_quadratic_max Day of maximal proximal treatment effect when beta_shape is quadratic.
#' @param tau_shape The pattern for expected availability; choices are constant, linear or quadratic.
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
#' @param tau_mean The average of expected availability.
#' @param tau_initial The initial Value of expected availability when tau_shape is linear or quadratic.
#' @param tau_quadratic_max The changing point of availability when tau_shape is quadratic.
#' @param dimB The number of parameters used in the main/average effect of proximal outcome.
#' @param power The desired power to achieve.
#' @param sigLev The significance level or type I error rate.
#'
#' @return The minimal sample size to achieve the desired power.
#'
#' @export calculateSampleSize
#' @references Seewald, N.J.; Sun, J.; Liao, P. "MRT-SS Calculator: An R Shiny Application for Sample Size
#' Calculation in Micro-Randomized Trials". arXiv:1609.00695
#' @examples
#' calculateSampleSize(days=42,
#'                     occ_per_day=5,
#'                     prob=0.4,
#'                     beta_shape="quadratic",
#'                     beta_mean=0.1,
#'                     beta_initial=0,
#'                     beta_quadratic_max=28,
#'                     tau_shape="quadratic",
#'                     tau_mean=0.5,
#'                     tau_initial=0.7,
#'                     tau_quadratic_max=42,
#'                     dimB=3,
#'                     power=0.8,
#'                     sigLev=0.05)
#'
#' prob1 <- c(replicate(35,0.7),replicate(35,0.6),replicate(35,0.5),replicate(35,0.4))
#' calculateSampleSize(days=28,
#'                     occ_per_day=5,
#'                     prob=prob1,
#'                     beta_shape="quadratic",
#'                     beta_mean=0.1,
#'                     beta_initial=0,
#'                     beta_quadratic_max=28,
#'                     tau_shape="quadratic",
#'                     tau_mean=0.5,
#'                     tau_initial=0.7,
#'                     tau_quadratic_max=42,
#'                     dimB=3,
#'                     power=0.8,
#'                     sigLev=0.05)

calculateSampleSize <- function(days,occ_per_day,prob,
                                beta_shape,beta_mean,beta_initial,beta_quadratic_max,
                                tau_shape,tau_mean,tau_initial,tau_quadratic_max,dimB,
                                power,sigLev){

  if(power <= 0){

    stop("Error: Please specify the power greater than 0")

  }
  if(power >= 1){

    stop("Error: Please specify the power less than 1")

  }
  beta_shape <- tolower(beta_shape)
  tau_shape <- tolower(tau_shape)
  validateParameters(days,occ_per_day,prob,
                     beta_shape,beta_mean,beta_initial,beta_quadratic_max,
                     tau_shape,tau_mean,tau_initial,tau_quadratic_max,dimB,
                     sigLev)

  beta_input <- generateBeta(days,occ_per_day,beta_shape,beta_mean,beta_initial,beta_quadratic_max)
  tau_input <- generateTau(days,occ_per_day,tau_shape,tau_mean,tau_initial,tau_quadratic_max)


  total <- days*occ_per_day
  input_avail <- vector('numeric', total)
  input_effect <- vector('numeric',total)

  p_input <- 3
  if(beta_shape == "constant"){
    p_input <- 1
  } else if(beta_shape == "linear"){
    p_input <- 2
  }

  for(k in 1:days){

    input_effect[(occ_per_day*k-occ_per_day+1):(occ_per_day*k)] = replicate(occ_per_day, beta_input[k])

  }


  for(k in 1:days){

    input_avail[(occ_per_day*k-occ_per_day+1):(occ_per_day*k)] = replicate(occ_per_day, tau_input[k])

  }

  if(length(prob) == 1 || length(prob) == days || length(prob) == days*occ_per_day){

      N <- SampleSize(days, occ_per_day, input_effect, input_avail, delta = prob, alpha0=sigLev, beta0=power, p=p_input, q=dimB, Nmax=1000)

      if(N > 10){

        print(sprintf("The required sample size is %3d to attain %5.2f%s power when the significance level is %4.2f.", N, power*100,"%",sigLev))

      }else{

        ### if the calculated sample size is less than 10, we won't output the exact sample size ###

        print(sprintf("The required sample size is less than or equal to 10 to attain  %5.2f %s power when the significance level is %4.2f.", power*100, "%", sigLev))

      }

  }else{

    print(paste("The input randomization probabilities (Prob) has wrong length!"))

  }
  return(N)

}


