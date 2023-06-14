rm(list = ls())
#load(file = "new-setup\\model1\\model1.RData")
set.seed(1111)
mu1_mean = 0; mu1_var = 0.0033
bet_mean = 7; bet_var = 3
q = seq(0.0001,0.9999, length.out = 40)
mu1 = qnorm(q, mu1_mean, sqrt(mu1_var)); 
#mu2 = qunif(q, -0.3, 0.3); 
bet = qnorm(q, bet_mean, sqrt(bet_var))

# cent = kmeans(cbind(mu1,bet),5)$centers
moment = c(mu1_mean, mu1_var + mu1_mean^2, 
           bet_mean, bet_var + bet_mean^2)

d = cbind(mu1 - moment[1],
          mu1^2 - moment[2],
          bet - moment[3],
          bet^2 - moment[4]
          )

fn = function(X,a) sum(exp(X%*%a)) 

lam = optim(c(0,0,0,0), fn = fn, X = d)$par
w = exp(d%*%lam)


Hmisc::wtd.mean(mu1,w); runif(10000,-0.1, 0.1)|> mean()
Hmisc::wtd.var(mu1,w); runif(10000,-0.1, 0.1)|> var()
Hmisc::wtd.mean(bet,w) ; runif(10000,4, 10)|> mean()
Hmisc::wtd.var(bet,w); runif(10000,4, 10)|> var()




wt = w
thet = cbind(mu = mu1, bet = bet)

#####
x = seq(-1,1,0.01)
p = rep(1/length(x), length(x))

emat = apply(thet,1, \(thet){
  
  lapply(x, \(xi){
    
    p = 1/(1+exp(-thet["bet"]*(xi-thet["mu"])))
    v = c(thet["bet"], thet["mu"] - xi)
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

save(m, file = "new-setup\\logistic2.RData") 

##################
################################
###############################################


rm(list = ls())
#load(file = "new-setup\\model1\\model1.RData")
set.seed(1111)
mu1_mean = 0; mu1_var = 0.03
bet_mean = 7; bet_var = 3
q = seq(0.0001,0.9999, length.out = 40)
mu1 = qnorm(q, mu1_mean, sqrt(mu1_var)); 
#mu2 = qunif(q, -0.3, 0.3); 
bet = qnorm(q, bet_mean, sqrt(bet_var))

# cent = kmeans(cbind(mu1,bet),5)$centers
moment = c(mu1_mean, mu1_var + mu1_mean^2, 
           bet_mean, bet_var + bet_mean^2)

d = cbind(mu1 - moment[1],
          mu1^2 - moment[2],
          bet - moment[3],
          bet^2 - moment[4])

fn = function(X,a) sum(exp(X%*%a)) 

lam = optim(c(0,0,0,0), fn = fn, X = d)$par
w = exp(d%*%lam)


Hmisc::wtd.mean(mu1,w); runif(10000,-0.3, 0.3)|> mean()
Hmisc::wtd.var(mu1,w); runif(10000,-0.3, 0.3)|> var()
Hmisc::wtd.mean(bet,w) ; runif(10000,4, 10)|> mean()
Hmisc::wtd.var(bet,w); runif(10000,4, 10)|> var()

wt = w
thet = cbind(mu = mu1, bet = bet)

#####
x = seq(-1,1,0.01)
p = rep(1/length(x), length(x))

emat = apply(thet,1, \(thet){
  
  lapply(x, \(xi){
    
    p = 1/(1+exp(-thet["bet"]*(xi-thet["mu"])))
    v = c(thet["bet"], thet["mu"] - xi)
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
save(m, file = "new-setup\\logistic3.RData") 

