library(fastGHS)
library(huge)
library(glasso)
library(igraph)
source('examples/GHS.r')

## Test fastGHS with AIC selection of taus


# EXAMPLE 1 ------------------------------------------------------------------
# GENERATE GRAPH with tau fixed: n=100, p=50, larger partial correlations (0.229)

n.fix=100
p.fix=50
set.seed(12345)
data.sf.fix = huge::huge.generator(n=n.fix, d=p.fix,graph = 'scale-free',v=0.5,u=0.05) 
g.true.sf.fix = data.sf.fix$theta # True adjacency matrix
theta.true.fix = data.sf.fix$omega # The precision matrix
theta.true.fix[which(theta.true.fix<10e-5,arr.ind=T)]=0  
g.sf.fix=graph.adjacency(data.sf.fix$theta,mode="undirected",diag=F) # true igraph object
x.sf.fix = data.sf.fix$data # Observed attributes. nxp matrix.
x.sf.scaled.fix= scale(x.sf.fix) # Scale columns/variables.
s.sf.scaled.fix = cov(x.sf.scaled.fix) # Empirical covariance matrix
data.sf.fix$sparsity # True sparsity: 0.04
# Look at precision matrix (partial correlations)
cov2cor(theta.true.fix[1:5,1:5])
#[,1]      [,2]      [,3]      [,4]      [,5]
#[1,] 1.0000000 0.2290578 0.0000000 0.0000000 0.0000000
#[2,] 0.2290578 1.0000000 0.2290578 0.0000000 0.0000000
#[3,] 0.0000000 0.2290578 1.0000000 0.2290578 0.2290578
#[4,] 0.0000000 0.0000000 0.2290578 1.0000000 0.0000000
#[5,] 0.0000000 0.0000000 0.2290578 0.0000000 1.0000000

system.time(
res.fix <- fastGHS(x.sf.fix,tau_sq = 0.1, AIC_selection = T, AIC_eps = 0.1, epsilon = 1e-3, fix_tau=TRUE)
)
theta.est.fix <- cov2cor(res.fix$theta)
theta.est.fix[which(abs(theta.est.fix) < 1e-5, arr.ind = T)] = 0
tailoredGlasso::sparsity(theta.est.fix!=0)
# 0.01632653
tailoredGlasso::precision(as.matrix(theta.true.fix!=0), theta.est.fix!=0)
# 1
tailoredGlasso::recall(as.matrix(theta.true.fix!=0), theta.est.fix!=0)
# 0.4081633
theta.est.fix[1:5,1:5]
#     [,1]      [,2]      [,3]      [,4]      [,5]
#[1,]    1 0.0000000 0.0000000 0.0000000 0.0000000
#[2,]    0 1.0000000 0.2495652 0.0000000 0.0000000
#[3,]    0 0.2495652 1.0000000 0.2188154 0.2536045
#[4,]    0 0.0000000 0.2188154 1.0000000 0.0000000
#[5,]    0 0.0000000 0.2536045 0.0000000 1.0000000

# Value of tau
res.fix$tau_sq
# 5.6001

# Compare to ordinary GHS, forced to same sparsity to allow for direct comparison
set.seed(123)
system.time(
  ghs.res.fix <- GHS(t(x.sf.fix)%*%x.sf.fix,n.fix,burnin=100,nmc=1000)
)
hist(ghs.res.fix$taus.samples) # Small tau, around size 1e-4
theta.est.ghs = cov2cor(apply(ghs.res.fix$thetas.sampled, c(1,2), mean))
theta.est.off.diag.ghs <- theta.est.ghs
diag(theta.est.off.diag.ghs) <- NA
# Look at the distribution of the precision matrix elements: no clear distinguishment between 
hist(c(theta.est.off.diag.ghs), breaks=100)
theta.est.ghs[which(abs(theta.est.ghs) < quantile(abs(theta.est.off.diag.ghs),1-tailoredGlasso::sparsity(theta.est.fix!=0),na.rm = T), arr.ind = T)] = 0
quantile(abs(theta.est.off.diag.ghs),1-tailoredGlasso::sparsity(theta.est.fix!=0),na.rm = T)
# 98.36735% 
# 0.1721756 
mean(ghs.res.fix$taus.samples)
# 0.0001879007
tailoredGlasso::sparsity(theta.est.ghs!=0)
#  0.01632653
tailoredGlasso::precision(theta.true.fix!=0,theta.est.ghs!=0)
# 1
tailoredGlasso::recall(theta.true.fix!=0,theta.est.ghs!=0)
# 0.4081633
theta.est.ghs[1:5,1:5]
#[,1]      [,2]      [,3]      [,4]      [,5]
#     [,1]      [,2]      [,3]      [,4]      [,5]
#[1,]    1 0.0000000 0.0000000 0.0000000 0.0000000
#[2,]    0 1.0000000 0.2281446 0.0000000 0.0000000
#[3,]    0 0.2281446 1.0000000 0.1833029 0.2419605
#[4,]    0 0.0000000 0.1833029 1.0000000 0.0000000
#[5,]    0 0.0000000 0.2419605 0.0000000 1.0000000

# We get the same precision and recall with the methods. 

# We also see that fastGHS is much faster the MCMC approach. 


# Look at their edge agreement
tailoredGlasso::confusion.matrix(theta.est.fix!=0,theta.est.ghs!=0)
#     [,1] [,2]
#[1,]   16    4
#[2,]    4 1201

# Compare to the graphical lasso by forcing it to the same sparsity as the ECM GHS esitmate
gg = huge(x.sf.fix,method='glasso', lambda=0.465)
gg$sparsity
# 0.01632653
tailoredGlasso::precision(as.matrix(theta.true.fix!=0), gg$icov[[1]]!=0)
# 0.8
tailoredGlasso::recall(as.matrix(theta.true.fix!=0), gg$icov[[1]]!=0)
# 0.3265306

# Better precision and recall for ECM GHS than the graphical lasso!



# EXAMPLE 2 ------------------------------------------------------------------
# GENERATE GRAPH with tau fixed: n=200, p=150, larger partial correlations (0.183)

n.fix2=200
p.fix2=150
set.seed(12345)
data.sf.fix2 = huge::huge.generator(n=n.fix2, d=p.fix2,graph = 'scale-free',v=0.5,u=0.05) 
g.true.sf.fix2 = data.sf.fix2$theta # True adjacency matrix
theta.true.fix2 = data.sf.fix2$omega # The precision matrix
theta.true.fix2[which(theta.true.fix2<10e-5,arr.ind=T)]=0  
g.sf.fix2=graph.adjacency(data.sf.fix2$theta,mode="undirected",diag=F) # true igraph object
x.sf.fix2 = data.sf.fix2$data # Observed attributes. nxp matrix.
x.sf.scaled.fix2= scale(x.sf.fix2) # Scale columns/variables.
s.sf.scaled.fix2 = cov(x.sf.scaled.fix2) # Empirical covariance matrix
data.sf.fix2$sparsity # True sparsity: 0.013333
# Look at precision matrix (partial correlations)
cov2cor(theta.true.fix2[1:5,1:5])
#[,1]      [,2]      [,3]      [,4]      [,5]
#[1,] 1.0000000 0.1826553 0.0000000 0.0000000 0.0000000
#[2,] 0.1826553 1.0000000 0.1826553 0.0000000 0.0000000
#[3,] 0.0000000 0.1826553 1.0000000 0.1826553 0.1826553
#[4,] 0.0000000 0.0000000 0.1826553 1.0000000 0.0000000
#[5,] 0.0000000 0.0000000 0.1826553 0.0000000 1.0000000

# Now the MCMC GHS fails, so we do not include any results from it here

res.fix2 <- fastGHS(x.sf.fix2,tau_sq = 0.1,AIC_selection = T,AIC_eps = 0.1, epsilon = 1e-3)
theta.est.fix2 <- cov2cor(res.fix2$theta)
theta.est.fix2[which(abs(theta.est.fix2) < 1e-5, arr.ind = T)] = 0
tailoredGlasso::sparsity(theta.est.fix2!=0)
# 0.00787472
tailoredGlasso::precision(as.matrix(theta.true.fix2!=0), theta.est.fix2!=0)
# 0.8295455
tailoredGlasso::recall(as.matrix(theta.true.fix2!=0), theta.est.fix2!=0)
# 0.4899329
theta.est.fix2[1:5,1:5]
#          [,1]      [,2]      [,3]     [,4]      [,5]
#[1,] 1.0000000 0.1953034 0.0000000 0.000000 0.0000000
#[2,] 0.1953034 1.0000000 0.2155343 0.000000 0.0000000
#[3,] 0.0000000 0.2155343 1.0000000 0.215094 0.1896857
#[4,] 0.0000000 0.0000000 0.2150940 1.000000 0.0000000
#[5,] 0.0000000 0.0000000 0.1896857 0.000000 1.0000000

# Value of tau
res.fix2$tau_sq
# 8.2001

# Compare to the graphical lasso by forcing it to the same sparsity as the ECM GHS esitmate
gg = huge(x.sf.fix2,method='glasso', lambda=0.268)
gg$sparsity
# 0.00787472
tailoredGlasso::precision(as.matrix(theta.true.fix2!=0), gg$icov[[1]]!=0)
# 0.6022727
tailoredGlasso::recall(as.matrix(theta.true.fix2!=0), gg$icov[[1]]!=0)
# 0.3557047

# Better precision and recall for ECM GHS than the graphical lasso!


# EXAMPLE 3 ------------------------------------------------------------------
# GENERATE GRAPH with tau fixed: n=300, p=400


n.fix3=300
p.fix3=400
set.seed(12345)
data.sf.fix3 = huge::huge.generator(n=n.fix3, d=p.fix3,graph = 'scale-free') 
g.true.sf.fix3 = data.sf.fix3$theta # True adjacency matrix
theta.true.fix3 = data.sf.fix3$omega # The precision matrix
theta.true.fix3[which(theta.true.fix3<10e-5,arr.ind=T)]=0  
g.sf.fix3=graph.adjacency(data.sf.fix3$theta,mode="undirected",diag=F) # true igraph object
x.sf.fix3 = data.sf.fix3$data # Observed attributes. nxp matrix.
x.sf.scaled.fix3= scale(x.sf.fix3) # Scale columns/variables.
s.sf.scaled.fix3 = cov(x.sf.scaled.fix3) # Empirical covariance matrix
data.sf.fix3$sparsity # True sparsity: 0.005
# Look at precision matrix (partial correlations)
cov2cor(theta.true.fix3[1:5,1:5])
#          [,1]      [,2]      [,3]      [,4]      [,5]
#[1,] 1.0000000 0.1530514 0.0000000 0.0000000 0.0000000
#[2,] 0.1530514 1.0000000 0.1530514 0.0000000 0.0000000
#[3,] 0.0000000 0.1530514 1.0000000 0.1530514 0.1530514
#[4,] 0.0000000 0.0000000 0.1530514 1.0000000 0.0000000
#[5,] 0.0000000 0.0000000 0.1530514 0.0000000 1.0000000

# Now the MCMC GHS fails, so we do not include any results from it here

res.fix3 <- fastGHS(x.sf.fix3,tau_sq = 0.01,AIC_selection = T, AIC_eps = 1, epsilon = 1e-3)
theta.est.fix3 <- cov2cor(res.fix3$theta)
theta.est.fix3[which(abs(theta.est.fix3) < 1e-5, arr.ind = T)] = 0
tailoredGlasso::sparsity(theta.est.fix3!=0)
# 0.01681704
tailoredGlasso::precision(as.matrix(theta.true.fix3!=0), theta.est.fix3!=0)
#  0.2041729
tailoredGlasso::recall(as.matrix(theta.true.fix3!=0), theta.est.fix3!=0)
# 0.6867168
theta.est.fix3[1:5,1:5]
#          [,1]      [,2]      [,3]      [,4]      [,5]
#[1,] 1.0000000 0.1769722 0.0000000 0.0000000 0.0000000
#[2,] 0.1769722 1.0000000 0.1348952 0.0000000 0.0000000
#[3,] 0.0000000 0.1348952 1.0000000 0.1636696 0.2119474
#[4,] 0.0000000 0.0000000 0.1636696 1.0000000 0.0000000
#[5,] 0.0000000 0.0000000 0.2119474 0.0000000 1.0000000

# Value of tau
res.fix3$tau_sq
# 2.8001

# Compare to the graphical lasso by forcing it to the same sparsity as the ECM GHS esitmate
gg = huge(x.sf.fix3,method='glasso', lambda=0.1464)
gg$sparsity
# 0.01681704
tailoredGlasso::precision(as.matrix(theta.true.fix3!=0), gg$icov[[1]]!=0)
# 0.2064083
tailoredGlasso::recall(as.matrix(theta.true.fix3!=0), gg$icov[[1]]!=0)
#  0.6942356

# Better precision and recall for ECM GHS than the graphical lasso!

# EXAMPLE 4 ------------------------------------------------------------------
# GENERATE GRAPH with tau fixed: n=100, p=200, smaller partial correlations (0.18)

n.fix4=100
p.fix4=200
set.seed(12345)
data.sf.fix4 = huge::huge.generator(n=n.fix4, d=p.fix4,graph = 'scale-free') 
g.true.sf.fix4 = data.sf.fix4$theta # True adjacency matrix
theta.true.fix4 = data.sf.fix4$omega # The precision matrix
theta.true.fix4[which(theta.true.fix4<10e-5,arr.ind=T)]=0  
g.sf.fix4=graph.adjacency(data.sf.fix4$theta,mode="undirected",diag=F) # true igraph object
x.sf.fix4 = data.sf.fix4$data # Observed attributes. nxp matrix.
x.sf.scaled.fix4= scale(x.sf.fix4) # Scale columns/variables.
s.sf.scaled.fix4 = cov(x.sf.scaled.fix4) # Empirical covariance matrix
data.sf.fix4$sparsity # True sparsity: 0.01
# Look at precision matrix (partial correlations)
cov2cor(theta.true.fix4[1:5,1:5])
#           [,1]      [,2]      [,3]      [,4]      [,5]
#[1,] 1.0000000 0.1807393 0.0000000 0.0000000 0.0000000
#[2,] 0.1807393 1.0000000 0.1807393 0.0000000 0.0000000
#[3,] 0.0000000 0.1807393 1.0000000 0.1807393 0.1807393
#[4,] 0.0000000 0.0000000 0.1807393 1.0000000 0.0000000
#[5,] 0.0000000 0.0000000 0.1807393 0.0000000 1.0000000

# Now the MCMC GHS fails, so we do not include any results from it here

res.fix4 <- fastGHS(x.sf.fix4,epsilon = 1e-3, AIC_selection = T, AIC_eps = 1)
theta.est.fix4 <- cov2cor(res.fix4$theta)
theta.est.fix4[which(abs(theta.est.fix4) < 1e-5, arr.ind = T)] = 0
tailoredGlasso::sparsity(theta.est.fix4!=0)
# 0.04512563
tailoredGlasso::precision(as.matrix(theta.true.fix4!=0), theta.est.fix4!=0)
# 0.1169265
tailoredGlasso::recall(as.matrix(theta.true.fix4!=0), theta.est.fix4!=0)
# 0.5276382
theta.est.fix4[1:5,1:5]
#[,1]      [,2]      [,3]      [,4]      [,5]
#[1,] 1.0000000 0.1173207 0.0000000 0.0000000 0.0000000
#[2,] 0.1173207 1.0000000 0.1803740 0.0000000 0.0000000
#[3,] 0.0000000 0.1803740 1.0000000 0.2611403 0.2457267
#[4,] 0.0000000 0.0000000 0.2611403 1.0000000 0.0000000
#[5,] 0.0000000 0.0000000 0.2457267 0.0000000 1.0000000

# Value of tau
res.fix4$tau_sq
# 3.2001

# Compare to the graphical lasso by forcing it to the same sparsity as the ECM GHS esitmate
gg = huge(x.sf.fix4,method='glasso', lambda=0.2879)
gg$sparsity
# 0.01005025
tailoredGlasso::precision(as.matrix(theta.true.fix4!=0), gg$icov[[1]]!=0)
# 0.35
tailoredGlasso::recall(as.matrix(theta.true.fix4!=0), gg$icov[[1]]!=0)
# 0.3517588

# Better precision and recall for ECM GHS than the graphical lasso!


# EXAMPLE 5 ------------------------------------------------------------------
# GENERATE GRAPH with tau fixed: n=50, p=100, larger partial correlations (0.192)

# Test sensitivity to initialisation

n.fix5=100
p.fix5=150
set.seed(12345)
data.sf.fix5 = huge::huge.generator(n=n.fix5, d=p.fix5,graph = 'scale-free',v=0.5,u=0.05) 
g.true.sf.fix5 = data.sf.fix5$theta # True adjacency matrix
theta.true.fix5 = data.sf.fix5$omega # The precision matrix
theta.true.fix5[which(theta.true.fix5<10e-5,arr.ind=T)]=0  
g.sf.fix5=huge::graph.adjacency(data.sf.fix5$theta,mode="undirected",diag=F) # true igraph object
x.sf.fix5 = data.sf.fix5$data # Observed attributes. nxp matrix.
x.sf.scaled.fix5= scale(x.sf.fix5) # Scale columns/variables.
s.sf.scaled.fix5 = cov(x.sf.scaled.fix5) # Empirical covariance matrix
data.sf.fix5$sparsity # True sparsity: 0.01333
# Look at precision matrix (partial correlations)
cov2cor(theta.true.fix5[1:5,1:5])
#[,1]      [,2]      [,3]      [,4]      [,5]
#[1,] 1.0000000 0.1826553 0.0000000 0.0000000 0.0000000
#[2,] 0.1826553 1.0000000 0.1826553 0.0000000 0.0000000
#[3,] 0.0000000 0.1826553 1.0000000 0.1826553 0.1826553
#[4,] 0.0000000 0.0000000 0.1826553 1.0000000 0.0000000
#[5,] 0.0000000 0.0000000 0.1826553 0.0000000 1.0000000

# First initialisation (none)

res.fix5 <- fastGHS(x.sf.fix5,tau_sq = 0.05,epsilon = 1e-3, AIC_selection = T, AIC_eps = 0.1)
theta.est.fix5 <- cov2cor(res.fix5$theta)
theta.est.fix5[which(abs(theta.est.fix5) < 1e-5, arr.ind = T)] = 0
tailoredGlasso::sparsity(theta.est.fix5!=0)
# 0.01243848
tailoredGlasso::precision(as.matrix(theta.true.fix5!=0), theta.est.fix5!=0)
# 0.3595122
tailoredGlasso::recall(as.matrix(theta.true.fix5!=0), theta.est.fix5!=0)
# 0.3355505
theta.est.fix5[1:5,1:5]
#[,1]      [,2]      [,3]      [,4]      [,5]
#[1,] 1.0000000  0.2388029 0.0000000 0.0000000  0.0000000
#[2,] 0.2388029  1.0000000 0.1585512 0.0000000 -0.2205855
#[3,] 0.0000000  0.1585512 1.0000000 0.2034902  0.0000000
#[4,] 0.0000000  0.0000000 0.2034902 1.0000000  0.0000000
#[5,] 0.0000000 -0.2205855 0.0000000 0.0000000  1.0000000

# Value of tau
res.fix5$tau_sq

# Second initialisation (random)

theta_init = matrix(runif(p.fix5^2,0,0.2),nrow=p.fix5)
diag(theta_init) = 1
theta_init = as.matrix(Matrix::nearPD(theta_init)$mat)
res.fix5.2 <- fastGHS(x.sf.fix5,theta=theta_init,tau_sq = 0.05,epsilon = 1e-3, AIC_selection = T, AIC_eps = 0.1)
theta.est.fix5.2 <- cov2cor(res.fix5.2$theta)
theta.est.fix5.2[which(abs(theta.est.fix5.2) < 1e-5, arr.ind = T)] = 0
tailoredGlasso::sparsity(theta.est.fix5.2!=0)
# 0.01259642
tailoredGlasso::precision(as.matrix(theta.true.fix5!=0), theta.est.fix5.2!=0)
# 0.3566434
tailoredGlasso::recall(as.matrix(theta.true.fix5!=0), theta.est.fix5.2!=0)
# 0.3422819
theta.est.fix5.2[1:5,1:5]
#[,1]      [,2]      [,3]      [,4]      [,5]
#[1,] 1.0000000  0.2388139 0.0000000 0.0000000  0.0000000
#[2,] 0.2388139  1.0000000 0.1586496 0.0000000 -0.2206028
#[3,] 0.0000000  0.1586496 1.0000000 0.2034593  0.0000000
#[4,] 0.0000000  0.0000000 0.2034593 1.0000000  0.0000000
#[5,] 0.0000000 -0.2206028 0.0000000 0.0000000  1.0000000

# Value of tau
res.fix5.2$tau_sq

# Very similar results as to those with no initialisation!


# Third initialisation (random)

theta_init = matrix(runif(p.fix5^2,0,0.1),nrow=p.fix5)
diag(theta_init) = 1
theta_init = as.matrix(Matrix::nearPD(theta_init)$mat)
res.fix5.3 <- fastGHS(x.sf.fix5,theta=theta_init,tau_sq = 0.05,epsilon = 1e-3, AIC_selection = T, AIC_eps = 1)
theta.est.fix5.3 <- cov2cor(res.fix5.3$theta)
theta.est.fix5.3[which(abs(theta.est.fix5.3) < 1e-5, arr.ind = T)] = 0
tailoredGlasso::sparsity(theta.est.fix5.3!=0)
# 0.01261545
tailoredGlasso::precision(as.matrix(theta.true.fix5!=0), theta.est.fix5.3!=0)
# 0.3685943
tailoredGlasso::recall(as.matrix(theta.true.fix5!=0), theta.est.fix5.3!=0)
# 0.3489933
theta.est.fix5.3[1:5,1:5]
#[,1]      [,2]      [,3]      [,4]      [,5]
#[1,]    1 0.0000000 0.0000000  0.0000000  0.0000000
#[2,]    0 1.0000000 0.2640438  0.0000000  0.0000000
#[3,]    0 0.2640438 1.0000000  0.2569356  0.0000000
#[4,]    0 0.0000000 0.2569356  1.0000000 -0.3453255
#[5,]    0 0.0000000 0.0000000 -0.3453255  1.0000000

# Value of tau
res.fix5.3$tau_sq

# Very similar results to the others, only a bit better. 
# Sparsity is preserved in all cases.


