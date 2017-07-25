#' @title  plot the graph for the proximal treatment effect
#'
#' @description plot of the graphs for the proximal treatment effect when the trend for
#' the proximal treatment effect is constant, linear or quadractic.
#'
#' @param days Duration of the study.
#' @param occ_per_day Number of decision time points per day.
#' @param beta_shape The trend for the proximal treatment effect, choices are constant, linear or quadratic.
#'  Note:
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
#'
#' @importFrom graphics plot points abline legend
#' @return A graph for the proximal treatment effect.
#'
#' @export plotProximalEffect
#' @examples
#' plotProximalEffect(days=42,
#'                    occ_per_day=5,
#'                    beta_shape="quadratic",
#'                    beta_mean=0.1,
#'                    beta_initial=0,
#'                    beta_quadratic_max=28)
#'
#'

plotProximalEffect <- function(days,occ_per_day,beta_shape,beta_mean,beta_initial,beta_quadratic_max){
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

  validateProximalEffectParameters(beta_shape,beta_mean,beta_initial)

  beta_shape <- tolower(beta_shape)
  beta_input <- generateBeta(days,occ_per_day,beta_shape,beta_mean,beta_initial,beta_quadratic_max)

  if(beta_shape == "constant"){
    up = 0.2
    if(max(beta_input) > 0.2)
    {
      up = max(beta_input)
    }
    plot(beta_input,xlab = "Days", ylab = "Proximal Effect", ylim = c(0,up), type = "o",
         pch = 16, cex = 0.8, col = 2)
    points(x = 1:length(beta_input), y = rep(0, length(beta_input)), type = 'o', pch = 16, col = 4, cex = 0.8)
    abline(h=beta_mean, lty = 2, col = 1)
    legend("topleft", cex = 1, legend=c('Null Hypothesis','Alternate Hypothesis','Average Effect'), col = c(4,2,1),lty = c(1,1,2), pch=c(16,16,NA),bty = "n")

  } else if(beta_shape == "linear"){
    if(min(beta_input) < 0){
      warning("Warning: Some values of the proximal treatment effect are less than 0.")
    }
    if(max(beta_input) > 1){
      warning("Warning: Some values of the proximal treatment effect are bigger than 1.")
    }

    up = 0.2
    if(max(beta_input) > 0.2)
    {
      up = max(beta_input)
    }
    plot(beta_input,xlab = "Days", ylab = "Proximal Effect", ylim = c(0,up), type = "o",
         pch = 16, cex = 0.8, col = 2)
    points(x = 1:length(beta_input), y = rep(0, length(beta_input)), type = 'o', pch = 16, col = 4, cex = 0.8)
    abline(h=beta_mean, lty = 2, col = 1)
    legend("topleft", cex = 1, legend=c('Null Hypothesis','Alternate Hypothesis','Average Effect'), col = c(4,2,1),lty = c(1,1,2), pch=c(16,16,NA),bty = "n")

  } else if(beta_shape == "quadratic"){
    if(min(beta_input) < 0){
      warning("Warning: Some values of the proximal treatment effect are less than 0.")
    }
    if(max(beta_input) > 1){
      warning("Warning: Some values of the proximal treatment effect are bigger than 1.")
    }
    up = 0.2
    if(max(beta_input) > 0.2)
    {
      up = max(beta_input)
    }

    plot(beta_input, xlab = "Days", ylab = "Effects", ylim = c(0,up), type = "o",
         pch = 16, cex = 0.8, col = 2)
    points(x = 1:length(beta_input), y = rep(0, length(beta_input)), type = 'o', pch = 16, col = 4, cex = 0.8)
    abline(h=beta_mean, lty = 2, col = 1)
    legend("topleft", cex = 1, legend=c('Null Hypothesis','Alternate Hypothesis','Average Effect'), col = c(4,2,1),lty = c(1,1,2), pch=c(16,16,NA),bty = "n")

  }

}

validateProximalEffectParameters <- function(beta_shape,beta_mean,beta_initial,beta_quadratic_max){
  if(beta_shape != "constant"
     && beta_shape !="linear"
     && beta_shape != "quadratic"){
    stop("Please select a pattern for proximal treatment effect")
  }

  if(beta_mean <= 0){
    stop("Error: Please specify an average standardized effect greater than 0")
  }

  if(beta_mean >= 1){
    stop("Error: Please specify an average standardized effect less than 1")
  }

  if(beta_shape == "quadratic"){
    if(beta_initial < 0){
      stop("Error: Please specify a Standardized Initial Effect greater than or equal to 0")
    }

    if(beta_initial <= beta_mean){
      stop("Error: Please specify the standardized Initial Effect less than or equal to average standardized effect")
    }
  }
}

getTrueBeta <- function(q1 = 0, q2 = 0.3, q3 = 28, days){

  ### output (coefficients of) standardized treatment effect or availability (quadratic
  ### form) given q1, q2, q3 ###
  ### only used when q3 (in days) is larger than half of the study

  beta = matrix(0, 3)
  beta[1] = q1
  a = sum(c(1:(days-1)))
  b = sum(c(1:(days-1))^2)
  mat = t(matrix(c(1, 2*(q3-1),a,b), 2, 2))
  beta[2:3] = solve(mat) %*% matrix(c(0, q2 * days))

  return(beta)
}


generateBeta <- function(days,occ_per_day,beta_shape,beta_mean,beta_initial,beta_quadratic_max){
  if(beta_shape == "quadratic"){
    ### Quadratic class of proximal treatment effect ###
    N <- days*occ_per_day
    H <- beta_quadratic_max * occ_per_day
    M <- beta_mean
    I <- beta_initial

    beta1 <- getTrueBeta(q1 = 0, q2 = M-I, q3 = beta_quadratic_max , days = days)
    b <- beta1[2]
    c <- beta1[3]
    sequence <- c(0:(days-1))
    beta_input <- I + b*sequence + c * sequence^2
  } else if(beta_shape == "constant"){

    ### Constant class of proximal treatment effect ###
    beta_input <- replicate(days, beta_mean)

  } else if(beta_shape == "linear"){

    ### Linear class of proximal treatment effect ###
    initial <- beta_initial
    mean <- beta_mean
    days <- days

    range <- (initial - mean) * 2

    if(days == 1) {
      beta_input <- mean
    } else if (range == 0) {
      beta_input <- replicate(days, mean)
    } else {
      num <- range / (days - 1)
      beta_input <- seq(from = mean + range / 2,
                        to = mean - range / 2,
                        by = -num)
    }
  }
  return(beta_input)
}

