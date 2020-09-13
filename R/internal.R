
#' @title  Internal MRT functions
#' @description   These internal functions are not usually called
#' directly by the user. Their names and capabilities may change
#' without warning from one version of \pkg{MRT} to the next.
#' @importFrom stats qf pf
#' @noRd
#'

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


generateTau <- function(days,occ_per_day,tau_shape,tau_mean,tau_initial,tau_quadratic_max){
  if(tau_shape == "quadratic"){
    ### Quadratic class of proximal treatment effect ###
    N <- days*occ_per_day
    H <- tau_quadratic_max * occ_per_day
    M <- tau_mean
    I <- 2*M - tau_initial

    beta1 <- getTrueBeta(q1 = 0, q2 = M-I, q3 = tau_quadratic_max , days = days)
    b <- beta1[2]
    c <- beta1[3]
    sequence = c(0:(days-1))
    tau_input <- I + b*sequence + c * sequence^2
    tau_input <- 2*M-tau_input
  } else if(tau_shape == "constant"){

    ### Constant class of expected availability ###
    tau_input <- replicate(days, tau_mean)

  } else if(tau_shape == "linear"){

    ### Linear class of expected availability ###
    initial <- tau_initial
    mean <- tau_mean
    days <- days

    range <- (initial - mean) * 2

    if(days == 1) {
      tau_input <- mean
    } else if (range == 0) {
      tau_input <- replicate(days, mean)
    } else {
      num <- range / (days - 1)
      tau_input <- seq(from = mean + range / 2,
                       to = mean - range / 2,
                       by = -num)
    }
  }
  return(tau_input)
}


validateParameters <- function(days,occ_per_day,prob,
                               beta_shape,beta_mean,beta_initial,beta_quadratic_max,
                               tau_shape,tau_mean,tau_initial,tau_quadratic_max,dimB,
                               sigLev){

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

  if(length(prob) == 1){
    if(prob <= 0){
      stop("Error: Please specify the randomization probability greater than 0")
    }
    if(prob >= 1){
      stop("Error: Please specify the randomization probability less than 1")
    }
  } else {
    if(length(prob) != days * occ_per_day)
    {
      stop("Error: The number of prob does not match the number of decision times provided in Study Setup.")
    }
    if(max(prob) > 1)
    {
      stop("Error: The provided randomization probability is greater than 1 for one or more decision times.")
    }
    if(min(prob) < 0)
    {
      stop("Error: The provided randomization probability is less than 0 for one or more decision times.")
    }
  }

  if(!((dimB == as.integer(dimB)) & (dimB > 0))){

    stop("Error: Please specify the parameters in the main effect to be a positive integer")

  }

  if(sigLev < 0){
    stop("Error: Please specify the significance level greater than or equal to 0")
  }
  if(sigLev > 1){
    stop("Error: Please specify the significance level less than or equal to 1")
  }

  validateProximalEffectParameters(beta_shape,beta_mean,beta_initial)

  validateEffectAvailParameters(tau_shape,tau_mean)


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

    if(beta_initial > beta_mean){
      stop("Error: Please specify the standardized Initial Effect less than or equal to average standardized effect")
    }
  }
}


validateEffectAvailParameters <- function(tau_shape,tau_mean){
  if(tau_shape != "constant"
     && tau_shape !="linear"
     && tau_shape != "quadratic"){
    stop("Please select a pattern for expected availability")
  }

  if(tau_mean <= 0){
    stop("Error: Please specify the mean of availability greatert than 0")
  }

}


para <- function(order, Total, days, occ_per_day){

  if(order == 0){

    return(matrix(1, Total, 1))

  }else{

    Z <- matrix(0, Total, order + 1)
    Z[,1] <- 1;
    for(i in 1:order){

      Z[, i+1] <- (rep(c(0:(days-1)), each =occ_per_day))^i

    }

  }

  return(Z)

}


PowerCal <- function(p, q, N, d, beta, alpha){

  C <- N * t(d) %*% beta %*% d ## noncentral para.
  df2 <- N-(p+q) ## degrees of freedom
  adj_c <- qf(1-alpha, df1 = p, df2 = df2) ### critical value

  # output the power
  power <- 1 - pf(q = adj_c, df1 = p, df2 = df2 , ncp = C)

  return(power)
}


SampleSize <- function(days, occ_per_day, beta_t, tau, delta, alpha0, beta0, p, q, Nmax=1000){

  ### MRT Sample Size Calculation (Liao et. al 2015) ###

  ### INPUT:
  ### standardized treatment effect (beta_t), expected availability (tau),
  ### randomization probability (delta), significance level (alpha0)
  ### desired power (beta0), study setup (days, occ.per.days)
  ### p, q = number of parameters in Z and B (assume Z are quadratic, linear or constant)
  ### Nmax is the maximum number of participants to acheive the power

  ### OUTPUT: Sample size


  Total <- days * occ_per_day;
  Z <- para(order = p-1, Total, days, occ_per_day);

  stopifnot(length(beta_t) == Total)
  stopifnot(length(tau) == Total);

  if(length(delta) == 1){

    delta <- rep(delta, Total);

  }else if(length(delta) == days){

    delta <- rep(delta, each = occ_per_day);

  }else{

    stopifnot(length(delta) == Total)
  }


  beta <- t(Z) %*% diag(tau * delta * (1-delta)) %*% Z

  d <- solve(t(Z) %*% diag(tau) %*% Z) %*% t(Z)  %*% diag(tau) %*% beta_t

  N.all <- c(10:Nmax);

  N <- Nmax
  for( k in N.all){
    power <- PowerCal(p,q,k,d,beta,alpha0)
    if(power >= beta0){
      N <- k
      break
    }
  }

  if(N == Nmax){

    stop(paste("Cannot attain",beta0,"power when sample size is below", Nmax))

  }

  return(N)
}


PowerCalculation <- function(days, occ_per_day, N, beta_t, tau, delta, alpha0, p, q){

  ### MRT Power Calculation (Liao et. al 2015) ###

  ### INPUT:
  ### sample size (N)
  ### standardized treatment effect (beta_t), expected availability (tau),
  ### randomization probability (delta), significance level (alpha0)
  ### p, q = number of parameters in Z and B (assume Z are quadratic, linear or constant)

  ### OUTPUT: (Estimated) power


  Total <- days * occ_per_day ;
  Z <- para(order = p-1, Total, days, occ_per_day);

  stopifnot(length(beta_t) == Total)
  stopifnot(length(tau) == Total);


  if(length(delta) == 1){

    delta <- rep(delta, Total);

  }else if(length(delta) == days){

    delta <- rep(delta, each = occ_per_day);

  }else{

    stopifnot(length(delta) == Total)
  }


  beta <- t(Z) %*% diag(tau * delta * (1-delta)) %*% Z

  d <- solve(t(Z) %*% diag(tau) %*% Z) %*% t(Z)  %*% diag(tau) %*% beta_t

  power <- PowerCal(p,q,N,d,beta,alpha0)

  return(power)
}

