CIRZeroBondPrice <- function(t, alpha, mu, sigma,lambda, x)
{
  # t es T-t
  gamma <- sqrt((alpha + lambda)^2 + 2 * sigma^2);
  d <- 2*gamma+(alpha+lambda+gamma)*(exp(t*gamma)-1)
  A <- 2*gamma*exp((alpha+lambda+gamma)*t/2)/d
  A <- A^(2* alpha*mu/sigma^2)
  B <- 2*(exp(t*gamma)-1)/d
  A*exp(-B*x)
}