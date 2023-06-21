
rm(list = ls())  
#load(file = "new-setup\\model1\\model1.RData")
set.seed(2222)
# beta ~ N(7,3), mu1 ~ N(0, 0.003), mu2 ~ N(0, 0.03)
beta_mean = 7; beta_var = 3
mu1_mean = 0; mu1_var = 0.003
#mu2_mean = 0; mu2_var = 0.03

l = 100
beta = qnorm(seq(0.0025,0.995, length.out = l), beta_mean, sqrt(beta_var)*2)
mu1 = qnorm(seq(0.0025,0.995, length.out = l), mu1_mean, sqrt(mu1_var)*2)
#mu2 = qnorm(seq(0.01,0.99, length.out = l), mu2_mean, sqrt(mu2_var))

d1 = cbind(beta - beta_mean, beta^2 - (beta_var + beta_mean^2),
           mu1 - mu1_mean, mu1^2 - (mu1_var +mu1_mean^2)
           )

fn = function(X,a) sum(exp(X%*%a)) 

lam = optim(c(0,0,0,0), fn = fn, X = d1)$par
w = exp(d1%*%lam)|> c()
wt = w/sum(w)

Hmisc::wtd.mean(beta,w) ;Hmisc::wtd.var(beta,w)
Hmisc::wtd.mean(mu1,w) ;Hmisc::wtd.var(mu1,w)

# 
a = -weighted.mean(log(wt),wt); a
# 
b =  (1/2*log(2*pi*beta_var) + 1/2 +  1/2*log(2*pi*mu1_var) + 1/2);b

a/b


f1 = rnorm(100000, beta_mean, sqrt(beta_var))|> 
  dnorm(beta_mean, sqrt(beta_var))
f2 = rnorm(100000, mu1_mean, sqrt(mu1_var))|> 
  dnorm(mu1_mean, sqrt(mu1_var))
-mean(log(f1*f2))




plot(0,0, ylim = c(0,0.25), xlim = c(min(beta), max(beta)))
segments(x0 = beta, y0 = 0, x1 = beta, y1 = w/sum(w), col = "red" )
points(beta,  w/sum(w), col = "red")

rnorm(100000, beta_mean, sqrt(beta_var))|> density()|> lines()






g = 100
curve(dnorm(x, mu1_mean, sqrt(mu1_var)), xlim = c(-0.3,0.3))
segments(x0 = mu1, y0 = 0, x1 = mu1, y1 = g*w/sum(w), col = "red" )
points(mu1,  g*w/sum(w), col = "red")
#-----------------------------------------

wt = w/sum(w)
thet = cbind(beta = beta, mu = mu1)


#####
x = seq(-1,1,0.01)
p = rep(1/length(x), length(x))

emat = apply(thet,1, \(thet){
  
  lapply(x, \(xi){
    
    p = 1/(1+exp(-thet["beta"]*(xi-thet["mu"])))
    v = c(thet["beta"], thet["mu"] - xi)
    v%*%t(v)*p*(1-p)
  })
}, simplify = FALSE)

iter = 0

repeat{
  iter = iter+1
  trM = 
    sapply(emat, \(m){
      M = mapply(`*`,m, p, SIMPLIFY = FALSE)|>
        {\(i) Reduce(`+`,i)}()
      sapply(m, \(i)  sum(diag(solve(M)%*%i)) )
    })
  
  d = matrixStats::rowWeightedMeans(trM,wt)
  
  p = (p*d^0.95)/sum(p*d^0.95)
  #plot(x,p, type = "o")
  #points(x,p, type = "o", col = "red")
  #print(max(d))
  #print(iter)
  if(iter == 50000) break
}
m = cbind(x, p, d)

save(m, file = "new-setup\\logistic-2-points.RData") 

##################
################################
###############################################




rm(list = ls())
#load(file = "new-setup\\model1\\model1.RData")
set.seed(1111)
# beta ~ N(7,3), mu1 ~ N(0, 0.003), mu2 ~ N(0, 0.03)
beta_mean = 7; beta_var = 3
mu1_mean = 0; mu1_var = 0.03
#mu2_mean = 0; mu2_var = 0.03

l = 50
beta = qnorm(seq(0.001,0.999, length.out = l), beta_mean, sqrt(beta_var))
mu1 = qnorm(seq(0.001,0.999, length.out = l), mu1_mean, sqrt(mu1_var))


#mu2 = qnorm(seq(0.01,0.99, length.out = l), mu2_mean, sqrt(mu2_var))

d1 = cbind(beta - beta_mean, beta^2 - (beta_var + beta_mean^2),
           beta^3 - (3*beta_mean*beta_var + beta_mean^3),
           mu1 - mu1_mean, mu1^2 - (mu1_var +mu1_mean^2),
           mu1^3 - (3*mu1_mean*mu1_var + mu1_mean^3)
           
)

fn = function(X,a) sum(exp(X%*%a)) 

lam = optim(c(0,0,0,0,0,0), fn = fn, X = d1)$par
w = exp(d1%*%lam)|> c()






Hmisc::wtd.mean(beta,w); Hmisc::wtd.var(beta,w); Weighted.Desc.Stat::w.skewness(beta,w)
rnorm(2000, beta_mean, sqrt(beta_var))|> 
  {\(i) c(mean(i), var(i),  e1071::skewness(i) ) }()


Hmisc::wtd.mean(mu1,w) ;Hmisc::wtd.var(mu1,w); Weighted.Desc.Stat::w.skewness(mu1,w)
rnorm(2000, mu1_mean, sqrt(mu1_var))|> 
  {\(i) c(mean(i), var(i),  e1071::skewness(i) ) }()

curve(dnorm(x, beta_mean, sqrt(beta_var*16)), xlim = c(-30,45) )
segments(x0 = beta, y0 = 0, x1 = beta, y1 = w/sum(w), col = "red" )
# points(beta,  w/sum(w), col = "red")



curve(dnorm(x, mu1_mean, sqrt(mu1_var*1555)) ,  xlim = c(-5,5))
segments(x0 = mu1, y0 = 0, x1 = mu1, y1 = w/sum(w), col = "red" )
#points(mu1,  w/sum(w), col = "red")
#-----------------------------------------

wt = w/sum(w)
thet = cbind(beta = beta, mu = mu1)


#####
x = seq(-1,1,0.01)
p = rep(1/length(x), length(x))

emat = apply(thet,1, \(thet){
  
  lapply(x, \(xi){
    
    p = 1/(1+exp(-thet["beta"]*(xi-thet["mu"])))
    v = c(thet["beta"], thet["mu"] - xi)
    v%*%t(v)*p*(1-p)
  })
}, simplify = FALSE)

iter = 0

repeat{
  iter = iter+1
  trM = 
    sapply(emat, \(m){
      M = mapply(`*`,m, p, SIMPLIFY = FALSE)|>
        {\(i) Reduce(`+`,i)}()
      sapply(m, \(i)  sum(diag(solve(M)%*%i)) )
    })
  
  d = matrixStats::rowWeightedMeans(trM,wt)
  
  p = (p*d^0.95)/sum(p*d^0.95)
  plot(x,p, type = "o")
  points(x,p, type = "o", col = "red")
  #print(max(d))
  #print(iter)
  if(iter == 50000) break
}
m = cbind(x, p, d)


plot(x,d, type = "o", col = "red")


save(m, file = "new-setup\\logistic-3-points.RData") 


##############
########################
########################################






rm(list = ls())
#load(file = "new-setup\\model1\\model1.RData")
set.seed(1111)
# beta ~ N(7,3), mu1 ~ N(0, 0.003), mu2 ~ N(0, 0.03)
beta_mean = 7; beta_var = 3
mu1_mean = 0; mu1_var = 1


beta = rnorm(1000, beta_mean, sqrt(beta_var))
mu1 = rnorm(1000, mu1_mean, sqrt(mu1_var))



#-----------------------------------------

thet = cbind(beta = beta, mu = mu1)
wt = rep(1/nrow(thet), nrow(thet) )


#####
x = seq(-1,1,0.01)
p = rep(1/length(x), length(x))

emat = apply(thet,1, \(thet){
  
  lapply(x, \(xi){
    
    p = 1/(1+exp(-thet["beta"]*(xi-thet["mu"])))
    v = c(thet["beta"], thet["mu"] - xi)
    v%*%t(v)*p*(1-p)
  })
}, simplify = FALSE)

iter = 0

repeat{
  iter = iter+1
  trM = 
    sapply(emat, \(m){
      M = mapply(`*`,m, p, SIMPLIFY = FALSE)|>
        {\(i) Reduce(`+`,i)}()
      sapply(m, \(i)  sum(diag(solve(M)%*%i)) )
    })
  
  d = matrixStats::rowWeightedMeans(trM,wt)
  
  p = (p*d^0.95)/sum(p*d^0.95)
  plot(x,p, type = "o")
  points(x,p, type = "o", col = "red")
  #print(max(d))
  #print(iter)
  if(iter == 50000) break
}
m = cbind(x, p, d)

save(m, file = "new-setup\\logistic-3-points_full.RData") 



