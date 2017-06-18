
#Newton Raphson or Fisher Scoring Method for Maximum Likelihood Estimation

x = c (1,9,1,5,6,8,2,4,2,8,7,7)
y = c (3,8,2,8,5,9,4,5,2,4,2,6)
N = length(y)
X = cbind (1,x)

#Initialize
theta <- c(2,3,4)
res0<- 1
to1 <- 1e-8
norm <- 10
iter <- 0
change <- 1

#Iterate
while (change>to1) {
grad.1 = (t(X) % * %y-t(X) % * % X% * %theta[1:2])/theta[3]
grad.2 = -N/ (2 * theta[3]) +1 / (2* theta[3] ^2) * (t(y-X% * % theta[1:2]) % * % (y - X% *%theta[1:2]))
g = rbind (grad.1, grad.2)
hess.1 = -(t(X)% * %X) / theta[3]
hess.2 = -N/(2 * theta[3]^2)
#Create Information matrix
p = ncol (X)
n = p+1
I = matrix(0,n,n)
I [1:p ,1:p] = hess.1
I[n*n] = hess.2
I= -I
theta = theta + solve(I)%*%g
Norm.1 = norm (theta, type = “F” )   #F = Frobenius; can also use type = “2”
Change <- abs(norm-norm.1)
norm = norm.1
iter <- iter+1
}
theta; iter; change; I

#Same Result Using R’s Optimization Function
loglike <- function(theta) { -sum(dnorm(y, mean = theta[1]+ theta[2] * x,  sd= sqrt(theta[3]), log = T))
}
maxmod <- optim (theta <- c(2,3,4), loglike, hessian = T, method = “BFGS”)
maxmod

