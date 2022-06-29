
fz <- function(r) log((1 + r) / (1 - r)) / 2 

fzi <- function(z){
  exp2z <- exp(2*z)
  (exp2z - 1) / (exp2z + 1) 
}

sez <- function(n) 1 / sqrt(n - 3)

zsc <- function(r, p0 = 0){
  (fz(r) - fz(p0)) * sqrt(n - 3)
}
