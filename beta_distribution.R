#http://eriqande.github.io/sisg_mcmc_course/s03-01-intro-mcmc-in-R.nb.html

likelihood = function(n,y,theta){
  return(theta^y*(1-theta)^(n-y))
}
theta = seq(from=0.01,to=0.99,by=0.01)

plot(theta, likelihood(400,72,theta))

loglike = function(n,y,theta){
  return(log(theta)*y+log(1-theta)*(n-y))
}
plot(theta, loglike(400,72,theta))



#problem
# number of solutions correct in an exam

theta = seq(from=0.00,to=1,by=0.01)
theta = seq(from=0.00,to=1,by=0.01)
plot(theta,dbeta(theta,1,1))
plot(theta,dbeta(theta,4,2)) #2/3rd
plot(theta,dbeta(theta,8,4)) #2/3

1- pbeta(0.25,8,4)
Student
plot(theta,dbeta(theta,8+33,4+40-33))


#p2
orig = read.table("http://www.randomservices.org/random/data/Challenger2.txt",header=T)



# Z = (observed - mean) / SD
#bionomial dist SD = sqrt(np(1-p)) b) mean = np

#Gamma theta ~ gamma(alpha, beta)


cs = read.table("~/Downloads/LM22.txt",sep="\t",header = T)
row.names(cs) = cs$Gene.symbol
cs = cs %>% dplyr::select(-Gene.symbol)
head(cs)
cs = cs[apply(cs,1,sd)>2500,]
library(pheatmap)
pheatmap(as.matrix(cs))


#### posterior update
prob = seq(0.1,0.9,0.1)
prior = c(0.06,0.06,0.06,0.06,0.52,0.06,0.06,0.06,0.06)
lhood = dbinom(4,20,prob = prob)
posterior  = prior * lhood / sum(prior * lhood)


## MCMC - Gamma distribution
m = 1000000 # no of sims
a = 2 #rshape 
b = 1/ 3 #rate 

theta = rgamma(m, a, b)
hist(theta)
curve(dgamma(x, a, b))

mean(theta)
a/b

#MC Standar error
se = sd(theta) / sqrt(m)

#conf -intervals 
mean(theta) - 2*se
mean(theta) + 2*se