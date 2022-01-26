library(fastGHS)
library(huge)
library(glasso)
library(igraph)
library(tailoredGlasso)
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
res.fix <- fastGHS(x.sf.fix, AIC_selection = T, AIC_eps = 0.1, epsilon = 1e-3)
)
theta.est.fix <- cov2cor(res.fix$theta)
theta.est.fix[which(abs(theta.est.fix) < 1e-5, arr.ind = T)] = 0
tailoredGlasso::sparsity(theta.est.fix!=0)
# 0.02040816
tailoredGlasso::precision(as.matrix(theta.true.fix!=0), theta.est.fix!=0)
# 1
tailoredGlasso::recall(as.matrix(theta.true.fix!=0), theta.est.fix!=0)
# 0.5102041
theta.est.fix[1:5,1:5]
#     [,1]      [,2]      [,3]      [,4]      [,5]
#[1,] 1.0000000 0.2187291 0.0000000 0.0000000 0.0000000
#[2,] 0.2187291 1.0000000 0.2796421 0.0000000 0.0000000
#[3,] 0.0000000 0.2796421 1.0000000 0.2360344 0.2375822
#[4,] 0.0000000 0.0000000 0.2360344 1.0000000 0.0000000
#[5,] 0.0000000 0.0000000 0.2375822 0.0000000 1.0000000

# Value of tau
res.fix$tau_sq
# 5.4001

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
# 97.95918% 
# 0.1407431
mean(ghs.res.fix$taus.samples)
# 0.000181281
tailoredGlasso::sparsity(theta.est.ghs!=0)
# 0.02040816
tailoredGlasso::precision(theta.true.fix!=0,theta.est.ghs!=0)
# 0.92
tailoredGlasso::recall(theta.true.fix!=0,theta.est.ghs!=0)
# 0.4693878
theta.est.ghs[1:5,1:5]
#[,1]      [,2]      [,3]      [,4]      [,5]
#[1,] 1.0000000 0.1707758 0.0000000 0.000000 0.0000000
#[2,] 0.1707758 1.0000000 0.2677585 0.000000 0.0000000
#[3,] 0.0000000 0.2677585 1.0000000 0.233041 0.2004593
#[4,] 0.0000000 0.0000000 0.2330410 1.000000 0.0000000
#[5,] 0.0000000 0.0000000 0.2004593 0.000000 1.0000000

# We get better precision and recall with ECMGHS. 

# We also see that fastGHS is much faster the MCMC approach. (over three times faster)


# Look at their edge agreement
tailoredGlasso::confusion.matrix(theta.est.fix!=0,theta.est.ghs!=0)
#     [,1] [,2]
# [1,]   23    2
# [2,]    2 1198

# Compare to the graphical lasso by forcing it to the same sparsity as the ECM GHS esitmate
gg = huge(x.sf.fix,method='glasso', lambda=0.46)
gg$sparsity
# 0.02040816
tailoredGlasso::precision(as.matrix(theta.true.fix!=0), gg$icov[[1]]!=0)
# 0.64
tailoredGlasso::recall(as.matrix(theta.true.fix!=0), gg$icov[[1]]!=0)
# 0.3265306

# Much better precision and recall for ECM GHS than the graphical lasso!



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

res.fix2 <- fastGHS(x.sf.fix2,AIC_selection = T,AIC_eps = 0.1, epsilon = 1e-3)
theta.est.fix2 <- cov2cor(res.fix2$theta)
theta.est.fix2[which(abs(theta.est.fix2) < 1e-5, arr.ind = T)] = 0
tailoredGlasso::sparsity(theta.est.fix2!=0)
# 0.007516779
tailoredGlasso::precision(as.matrix(theta.true.fix2!=0), theta.est.fix2!=0)
# 0.8690476
tailoredGlasso::recall(as.matrix(theta.true.fix2!=0), theta.est.fix2!=0)
# 0.4899329
theta.est.fix2[1:5,1:5]
#          [,1]      [,2]      [,3]     [,4]      [,5]
#[1,] 1.0000000 0.1771121 0.0000000 0.000000 0.0000000
#[2,] 0.1771121 1.0000000 0.2108469 0.000000 0.0000000
#[3,] 0.0000000 0.2108469 1.0000000 0.212593 0.2175741
#[4,] 0.0000000 0.0000000 0.2125930 1.000000 0.0000000
#[5,] 0.0000000 0.0000000 0.2175741 0.000000 1.0000000

# Value of tau
res.fix2$tau_sq
# 8.2001

# Compare to the graphical lasso by forcing it to the same sparsity as the ECM GHS esitmate
gg = huge(x.sf.fix2,method='glasso', lambda=0.2785)
gg$sparsity
# 0.007516779
tailoredGlasso::precision(as.matrix(theta.true.fix2!=0), gg$icov[[1]]!=0)
# 0.6190476
tailoredGlasso::recall(as.matrix(theta.true.fix2!=0), gg$icov[[1]]!=0)
# 0.3489933

# Much better precision and recall for ECM GHS than the graphical lasso!


# EXAMPLE 3 ------------------------------------------------------------------
# GENERATE GRAPH with tau fixed: n=300, p=400, medium partial cors (0.158)

n.fix3=300
p.fix3=400
set.seed(12345)
data.sf.fix3 = huge::huge.generator(n=n.fix3, d=p.fix3,graph = 'scale-free', v=1,u=0.0001) 
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
#[1,] 1.0000000 0.1578818 0.0000000 0.0000000 0.0000000
#[2,] 0.1578818 1.0000000 0.1578818 0.0000000 0.0000000
#[3,] 0.0000000 0.1578818 1.0000000 0.1578818 0.1578818
#[4,] 0.0000000 0.0000000 0.1578818 1.0000000 0.0000000
#[5,] 0.0000000 0.0000000 0.1578818 0.0000000 1.0000000

# Now the MCMC GHS fails, so we do not include any results from it here. 
# Error we get from GHS:
#Error in mu_i + solve(inv_C_chol, rnorm(p - 1)) : 
#  non-numeric argument to binary operator
#Timing stopped at: 3323 91.71 3414

res.fix3 <- fastGHS(x.sf.fix3,AIC_selection = T, AIC_eps = 1, epsilon = 1e-3)
theta.est.fix3 <- cov2cor(res.fix3$theta)
theta.est.fix3[which(abs(theta.est.fix3) < 1e-5, arr.ind = T)] = 0
tailoredGlasso::sparsity(theta.est.fix3!=0)
# 0.01612782
tailoredGlasso::precision(as.matrix(theta.true.fix3!=0), theta.est.fix3!=0)
# 0.2253302
tailoredGlasso::recall(as.matrix(theta.true.fix3!=0), theta.est.fix3!=0)
# 0.726817
theta.est.fix3[1:5,1:5]
#          [,1]      [,2]      [,3]      [,4]      [,5]
#[1,] 1.0000000 0.1167086 0.0000000 0.0000000 0.0000000
#[2,] 0.1167086 1.0000000 0.1518077 0.0000000 0.0000000
#[3,] 0.0000000 0.1518077 1.0000000 0.1715039 0.1187755
#[4,] 0.0000000 0.0000000 0.1715039 1.0000000 0.0000000
#[5,] 0.0000000 0.0000000 0.1187755 0.0000000 1.0000000

# Value of tau
res.fix3$tau_sq
# 6.4001

# Compare to the graphical lasso by forcing it to the same sparsity as the ECM GHS esitmate
gg = huge(x.sf.fix3,method='glasso', lambda=0.15775)
gg$sparsity
# 0.01612782
tailoredGlasso::precision(as.matrix(theta.true.fix3!=0), gg$icov[[1]]!=0)
# 0.2144522
tailoredGlasso::recall(as.matrix(theta.true.fix3!=0), gg$icov[[1]]!=0)
#  0.6917293

# Better precision and recall for ECM GHS than the graphical lasso!

# EXAMPLE 4 ------------------------------------------------------------------
# GENERATE GRAPH with tau fixed: n=100, p=200, medium partial correlations (0.18)

n.fix4=100
p.fix4=200
set.seed(12345)
data.sf.fix4 = huge::huge.generator(n=n.fix4, d=p.fix4,graph = 'scale-free', v=10, u=0.001) 
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
#         [,1]     [,2]     [,3]     [,4]     [,5]
#[1,] 1.000000 0.190733 0.000000 0.000000 0.000000
#[2,] 0.190733 1.000000 0.190733 0.000000 0.000000
#[3,] 0.000000 0.190733 1.000000 0.190733 0.190733
#[4,] 0.000000 0.000000 0.190733 1.000000 0.000000
#[5,] 0.000000 0.000000 0.190733 0.000000 1.000000

# Now the MCMC GHS fails, so we do not include any results from it here

# Use larger epsilon for AIC due to the large number of variables

res.fix4 <- fastGHS(x.sf.fix4,epsilon = 1e-5, AIC_selection = T, AIC_eps = 2)
theta.est.fix4 <- cov2cor(res.fix4$theta)
theta.est.fix4[which(abs(theta.est.fix4) < 1e-5, arr.ind = T)] = 0
tailoredGlasso::sparsity(theta.est.fix4!=0)
# 0.04623116
tailoredGlasso::precision(as.matrix(theta.true.fix4!=0), theta.est.fix4!=0)
# 0.1233596
tailoredGlasso::recall(as.matrix(theta.true.fix4!=0), theta.est.fix4!=0)
# 0.4723618
theta.est.fix4[1:5,1:5]
#[,1]      [,2]      [,3]      [,4]      [,5]
#[1,] 1.0000000 0.1696725 0.0000000 0.0000000 0.00000
#[2,] 0.1696725 1.0000000 0.3222101 0.0000000 0.00000
#[3,] 0.0000000 0.3222101 1.0000000 0.3955197 0.34098
#[4,] 0.0000000 0.0000000 0.3955197 1.0000000 0.00000
#[5,] 0.0000000 0.0000000 0.3409800 0.0000000 1.00000

# Value of tau
res.fix4$tau_sq
# 4.2001

# Compare to the graphical lasso by forcing it to the same sparsity as the ECM GHS esitmate
gg = huge(x.sf.fix4,method='glasso', lambda=0.5688)
gg$sparsity
# 0.03829146
tailoredGlasso::precision(as.matrix(theta.true.fix4!=0), gg$icov[[1]]!=0)
# 0.07874016
tailoredGlasso::recall(as.matrix(theta.true.fix4!=0), gg$icov[[1]]!=0)
# 0.3015075

# Better precision and recall for ECM GHS than the graphical lasso!


# EXAMPLE 5 ------------------------------------------------------------------
# GENERATE GRAPH with tau fixed: n=100, p=150

# Test sensitivity to initialisation

n.fix5=100
p.fix5=150
set.seed(12345)
data.sf.fix5 = huge::huge.generator(n=n.fix5, d=p.fix5,graph = 'scale-free',v=0.5,u=0.05) 
g.true.sf.fix5 = data.sf.fix5$theta # True adjacency matrix
theta.true.fix5 = data.sf.fix5$omega # The precision matrix
theta.true.fix5[which(theta.true.fix5<10e-5,arr.ind=T)]=0  
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

res.fix5 <- fastGHS(x.sf.fix5,epsilon = 1e-3, AIC_selection = T, AIC_eps = 0.1)
theta.est.fix5 <- cov2cor(res.fix5$theta)
theta.est.fix5[which(abs(theta.est.fix5) < 1e-5, arr.ind = T)] = 0
tailoredGlasso::sparsity(theta.est.fix5!=0)
# 0.02478747
tailoredGlasso::precision(as.matrix(theta.true.fix5!=0), theta.est.fix5!=0)
# 0.234657
tailoredGlasso::recall(as.matrix(theta.true.fix5!=0), theta.est.fix5!=0)
# 0.4362416
theta.est.fix5[1:5,1:5]
#[,1]      [,2]      [,3]      [,4]      [,5]
#[1,] 1.0000000  0.1581617 0.0000000 0.000000  0.0000000
#[2,] 0.1581617  1.0000000 0.1715320 0.000000 -0.1875500
#[3,] 0.0000000  0.1715320 1.0000000 0.207298  0.1284808
#[4,] 0.0000000  0.0000000 0.2072980 1.000000  0.0000000
#[5,] 0.0000000 -0.1875500 0.1284808 0.000000  1.0000000

# Value of tau
res.fix5$tau_sq
# 2.8001

# Second initialisation (random)

theta_init = matrix(runif(p.fix5^2,0,0.2),nrow=p.fix5)
diag(theta_init) = 1
theta_init = as.matrix(Matrix::nearPD(theta_init)$mat)
res.fix5.2 <- fastGHS(x.sf.fix5,theta=theta_init, epsilon = 1e-3, AIC_selection = T, AIC_eps = 0.1)
theta.est.fix5.2 <- cov2cor(res.fix5.2$theta)
theta.est.fix5.2[which(abs(theta.est.fix5.2) < 1e-5, arr.ind = T)] = 0
tailoredGlasso::sparsity(theta.est.fix5.2!=0)
# 0.02478747
tailoredGlasso::precision(as.matrix(theta.true.fix5!=0), theta.est.fix5.2!=0)
# 0.234657
tailoredGlasso::recall(as.matrix(theta.true.fix5!=0), theta.est.fix5.2!=0)
# 0.4362416
theta.est.fix5.2[1:5,1:5]
#          [,1]       [,2]      [,3]      [,4]       [,5]
#[1,] 1.0000000  0.1581615 0.0000000 0.0000000  0.0000000
#[2,] 0.1581615  1.0000000 0.1715322 0.0000000 -0.1875501
#[3,] 0.0000000  0.1715322 1.0000000 0.2072975  0.1284808
#[4,] 0.0000000  0.0000000 0.2072975 1.0000000  0.0000000
#[5,] 0.0000000 -0.1875501 0.1284808 0.0000000  1.0000000

# Value of tau
res.fix5.2$tau_sq
# 2.8001

# Same results as those with no initialisation!


# Third initialisation (random)m smaller values

theta_init = matrix(runif(p.fix5^2,0,0.1),nrow=p.fix5)
diag(theta_init) = 1
theta_init = as.matrix(Matrix::nearPD(theta_init)$mat)
res.fix5.3 <- fastGHS(x.sf.fix5,theta=theta_init, epsilon = 1e-3, AIC_selection = T, AIC_eps = 0.1)
theta.est.fix5.3 <- cov2cor(res.fix5.3$theta)
theta.est.fix5.3[which(abs(theta.est.fix5.3) < 1e-5, arr.ind = T)] = 0
tailoredGlasso::sparsity(theta.est.fix5.3!=0)
# 0.02478747
tailoredGlasso::precision(as.matrix(theta.true.fix5!=0), theta.est.fix5.3!=0)
# 0.234657
tailoredGlasso::recall(as.matrix(theta.true.fix5!=0), theta.est.fix5.3!=0)
# 0.4362416
theta.est.fix5.3[1:5,1:5]
#[,1]      [,2]      [,3]      [,4]      [,5]
#[1,] 1.0000000  0.1581617 0.0000000 0.000000  0.0000000
#[2,] 0.1581617  1.0000000 0.1715321 0.000000 -0.1875500
#[3,] 0.0000000  0.1715321 1.0000000 0.207298  0.1284808
#[4,] 0.0000000  0.0000000 0.2072980 1.000000  0.0000000
#[5,] 0.0000000 -0.1875500 0.1284808 0.000000  1.0000000

# Value of tau
res.fix5.3$tau_sq
# 2.8001

# Same results as the others.
# Initialisation did not matter here.

# Fourth initialisation (truth)

res.fix5.4 <- fastGHS(x.sf.fix5,theta=as.matrix(round(theta.true.fix5, 8)), epsilon = 1e-5, AIC_selection = T, AIC_eps = 0.1)
theta.est.fix5.4 <- cov2cor(res.fix5.4$theta)
theta.est.fix5.4[which(abs(theta.est.fix5.4) < 1e-5, arr.ind = T)] = 0
tailoredGlasso::sparsity(theta.est.fix5.4!=0)
# 0.02478747
tailoredGlasso::precision(as.matrix(theta.true.fix5!=0), theta.est.fix5.4!=0)
# 0.234657
tailoredGlasso::recall(as.matrix(theta.true.fix5!=0), theta.est.fix5.4!=0)
# 0.4362416
theta.est.fix5.4[1:5,1:5]
#[,1]      [,2]      [,3]      [,4]      [,5]
#[1,] 1.0000000  0.1581617 0.0000000 0.000000  0.0000000
#[2,] 0.1581617  1.0000000 0.1715321 0.000000 -0.1875500
#[3,] 0.0000000  0.1715321 1.0000000 0.207298  0.1284808
#[4,] 0.0000000  0.0000000 0.2072980 1.000000  0.0000000
#[5,] 0.0000000 -0.1875500 0.1284808 0.000000  1.0000000

# Value of tau
res.fix5.4$tau_sq
# 2.8001
