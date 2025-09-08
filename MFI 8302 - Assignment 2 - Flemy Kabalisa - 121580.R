#-------------------------- MFI 8302 - Assignment 2 ----------------------------------#

#---------------------- Flemy Kabalisa - 121580 ---------------------------#


#--------------------------------------------#

#------- Question 1 - Monte Carlo pricing of a European call ----------#

#--------------------------------------------#


#-----------------#

#--- Part a ----#

#-----------------#


#We set up our parameters

s0 <- 100
strike1 <- 100
r <- 0.05
sigma <- 0.2
T1 <- 0.5

#we discretise our parameters

Nsteps <- 252
Ms <- c(100,500,1000,5000,10000)

set.seed(12345)

#we set up the Black Scholes formula for comparison

bs_call <- function(s0, strike1, r, sigma, T1) {
  if (T1<= 0) return(max(s0-strike1, 0))
  d1 <- (log(s0/strike1) + (r + 0.5 *sigma^2) *T1) / (sigma* sqrt(T1))
  d2 <- d1 - sigma * sqrt(T1)
  price <- s0 * pnorm(d1) - strike1 * exp(-r * T1) * pnorm(d2)
  return(price)
}

#Euler simulation for S_T

simulate_ST_euler <- function(M, Nsteps, s0, r, sigma, T1) {
  dt <- T1/Nsteps
  S <- rep(s0,M)
  #we only need S_T, not full paths
  for (k in seq_len(Nsteps)){
    Z <- rnorm(M)
    S <- S + r * S * dt + sigma * S * sqrt(dt) * Z
  }
  return(S)
}

#Exact sampling from closed-form GBM

simulate_ST_exact <- function(M, s0, r, sigma, T1){
  Z <- rnorm(M)
  ST <- s0 * exp((r - 0.5 * sigma^2) * T1 + sigma * sqrt(T1) * Z)
  return(ST)
}

#Monte Carlo estimator for discounted payoffs

mc_price_from_ST <- function(ST, strike1, r, T1){
  payoff <- pmax(ST - strike1, 0)
  price <- exp(-r * T1) * mean(payoff)
  se <- exp(-r * T1) * sd(payoff) / sqrt(length(payoff))
  ci95 <- price + c(-1.96, 1.96) * se
  return(list(price = price, se = se, ci95 = ci95))
}

#--- Run simulations for the different M values 

res_euler <- data.frame(M = Ms, price = NA, se = NA, l95 = NA, u95 = NA)
res_exact <- data.frame(M = Ms, price = NA, se = NA, l95 = NA, u95 = NA)

for (i in seq_along(Ms)){
  M <- Ms[i]
  
  #Euler discretisation
  ST_e <- simulate_ST_euler(M, Nsteps, s0, r, sigma, T1)
  out_e <- mc_price_from_ST(ST_e, strike1, r, T1)
  res_euler[i, c("price", "se", "l95", "u95")] <- c(out_e$price, out_e$se, out_e$ci95)
  
  #Exact sampling
  ST_x <- simulate_ST_exact(M, s0, r, sigma, T1)
  out_x <- mc_price_from_ST(ST_x, strike1, r, T1)
  res_exact[i, c("price", "se", "l95", "u95")] <- c(out_x$price, out_x$se, out_x$ci95)
}

#True Black-Scholes price for comparison

bs_true <- bs_call(s0, strike1, r, sigma, T1)

#Summary Tables
cat("\nBlack Scholes (true) call price:", round(bs_true, 6), "\n\n")

cat("Euler discretisation results:\n")
print(res_euler)

cat("\nExact-sampling results:\n")
print(res_exact)


#--- Now we creat a plot to compare the convergence ---#

ylim_low <- min(res_euler$l95, res_exact$l95, bs_true) *0.98
ylim_up <- max(res_euler$u95, res_exact$u95, bs_true) *1.02

plot(Ms, res_euler$price, type = 'b', pch= 16, lty=1, log='x',
     xlab = "Number of Monte Carlo paths M (log scale)", ylab = "Estimated call price",
     ylim = c(ylim_low, ylim_up), main = "MC prices (Euler vs exact sampling) and BS true price")
arrows(Ms, res_euler$l95, Ms, res_euler$u95, angle = 90, code=3, length = 0.03)

points(Ms, res_exact$price, type = 'b', pch = 17, lty = 2, col="blue")
arrows(Ms, res_exact$l95, Ms, res_exact$u95, angle = 90, code = 3, length = 0.03, col = "blue")

abline(h = bs_true, col= "red", lwd = 2) #Our analytic BS price
legend("bottomright", legend = c("Euler MC", "Exact-sampling MC", "Black-Scholes (true)"),
       pch= c(16,17,NA), lty = c(1,2,1), col = c("black", "blue", "red"), bty="n")
