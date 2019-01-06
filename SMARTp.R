#### Sample Size calculation for the "SMARTp" Design from "SMARTp: A SMART design for non-surgical treatments of chronic periodontitis with spatially-referenced and non-randomly 
#### missing skewed outcomes", by Xu, Bandyopadhyay, Mirzaei, Michalowicz, and Chakraborty, (Under Review)


#### A Summary of the Design 

## SMART design with two stages, and 8 possible treatments (please see draft manuscript)
## Stage-1 include two treatments, i.e. treatments 3 and 8
## Patients who respond to the stage-1 treatment will receive same treatment at stage-2, while non-responders will be allocated to other treatments, i.e. randomly 
## allocated to treatments 4-7 at stage-2. Hence, there are 8 regimes for this design. 
## Outcome measures are continuous, clustered, and skewed, modeled via a skew-t density 
## The covariance structure within a cluster is spatially-referenced, modeled via a conditionally autoregressive (CAR) structure
## Each cluster sub-unit has a binary missingness indicator, which is associated to the corresponding continuous outcome measure through a shared-parameter (joint) model  



#### Sample size calculation is implemented via the function "SampleSize.SMARTp()" (see below) for 

## (a) Detecting effect size of a single treatment regime, and/or 
## (b) Difference between two treatment regimes with or without sharing initial treatment



#### We first list some necessary functions

#The CAR covariance structure 


CAR.Cov.Teeth=function(m, rho, tau)
{
M_m=matrix(0,m,m)
D_m=matrix(0,m,m)
if(m>1)
{
for(i in (1:(m-1)))
{
D_m[i, i+1]=1;D_m[i+1,i]=1;
}
diag(D_m)=1
diag(M_m)=3;M_m[1,1]=2;M_m[m,m]=2
Sigma=tau^2*solve(M_m-rho*D_m)
}
else{Sigma=tau^2}
return(Sigma)
}

#Estimated variance of \bar{Y_{i}}=\Sum_{t=1}^{28} Y_{it} M_{it} / \Sum_{t=1}^{28} M_{it} by Monte Carlo method

MC.Var.Yibar.Mis=function(mu, Sigma, sigma1, lambda, nu, sigma0, Num, a0, b0, cutoff)
{
m=dim(Sigma)[1]
Qit=rmvnorm(Num,matrix(0,1,m),as.matrix(Sigma))
Yit=matrix(1,Num,1)%*%matrix(mu, 1,m) + Qit + matrix( rst( n=Num*m, xi=0, omega=sigma1, alpha=lambda, nu=nu ),Num,m )
Iit=a0 + Qit*b0 + rmvnorm( Num,matrix(0,1,m),as.matrix( diag( rep( sigma0^2,m ) ) ) )
Mit=ifelse( Iit>cutoff,1,0 )
mYi=mean( rowSums( Yit*( 1-Mit ) )/rowSums( 1-Mit ),na.rm=T )
VarYi=var( rowSums( Yit*( 1-Mit ) )/rowSums( 1-Mit ),na.rm=T )
PM=mean(rowSums( Mit )/m)
#Note that available teeth have higher Yit than missing teeth.
return(list(Yit=Yit, Mit=Mit, Iit=Iit, mYi=mYi, VarYi=VarYi, PM=PM))
}



b0fun=function(c_i, b0, Sigma, sigma1, nu, lambda, sigma0) 
{
if(nu<Inf)
{
mean( b0*diag(Sigma)/sqrt( ( diag(Sigma) + ( sigma1^2*nu/( nu-2 ) - nu/pi*( gamma(0.5*(nu-1))/gamma(0.5*nu) )^2*sigma1^2*( lambda^2/( 1+lambda^2 ) ) ) )*( b0^2 * diag(Sigma) + sigma0^2 ) ) )-c_i 
}
else
{
mean( b0*diag(Sigma)/sqrt( ( diag(Sigma) + ( sigma1^2 - 2/pi*sigma1^2*( lambda^2/( 1+lambda^2 ) ) ) )*( b0^2 * diag(Sigma) + sigma0^2 ) ) )-c_i 
}
}


pifun=function(cutoff, a0, b0, Sigma, sigma0)
{
Epit=rep(0, 28)
for(j in 1:28)
{
Epit[j]=pnorm( (cutoff-a0)/sqrt( b0^2 * Sigma[j,j] + sigma0^2 ) )
}
mean(Epit)
}

a0fun=function(p_i, cutoff, a0, b0, Sigma, sigma0)
{
pifun(cutoff, a0, b0, Sigma, sigma0)-p_i
}



### SAMPLE SIZE CALCULATION FUNCTION


SampleSize.SMARTp = function(mu, st1, dtr, regime, pow, b, a, rho, tau, sigma1, lambda, nu, sigma0, Num, p_i, c_i, a0, b0, cutoff) 
{
#mu: mean matrix of dimension (# of treatment path X m), where a row represents a corresponding treatment path, a column represents a corresponding unit within a cluster

#st1: stage-1 treatment matrix, where rows represent the corresponding stage-1 treatments, the first column includes the numbers of treatment options for responder, 
#     the second column includes the numbers of treatment options for non-responder, the third column are the responding rates, and the fourth column includes the row numbers

#dtr: no.of dynamic treatment regime (DTR) by four matrix, the first column represents the number of DTRs, the second column represents the correpsonding number of treatment path if responds for the correponding DTRs in the first column, the third column represents the corresponding number of treatment path if non-responds for the corresponding DTRs in the first column, the fourth column represents the corresponding initial treatment

#regime: treatment regime scalar or vector, i.e. length should be 1 or 2



#pow=power 
#b=beta, type II error rate
#a=alpha, type I error rate,


#tau:    variation paramter of the CAR model
#rho:    association parameter of the CAR model
#sigma1: standard deviation of the residual for the continuous Y_{it}
#lambda: skewness parameter of the residual for the continuous Y_{it}
#nu:     degrees of freedom, or kurtosis parameter from the residual for the continuous outcome Y_{it}

# Note, the residual for the continuous outcome Y_{it} can follow normal, skew-normal, t or skew-t distributions
# sigma0: standard deviation of the residual for the binary outcome M_{it}
# Num:    number of samples to estimate variance of \bar{Y}_i


#p_i: the expected proportion of available teeth for patient 'i', i.e p_i=E(\sum_{t=1}^{28}(1-M_{it})/28)  
#c_i: the average Pearson correlation coefficient between Y_{it} and M_{it0} over the 28 teeth, i.e. t=1,...,28, for patient 'i'
#a0: intercept parameter in the probit regression model for the binary outcome M_{it}
#b0: slope parameter corresponding to the spatial random effect in the probit regression model for the binary outcome M_{it}
#cutoff: cut-off value of the binary outcome regression


if(missing(pow)){pow=0.8} 
if(missing(b)){b=0.2}
if(missing(a)){a=0.05}
if(missing(tau)){tau=0.85}
if(missing(rho)){rho=0.975}
if(missing(sigma1)){sigma1=0.95}
if(missing(sigma0)){sigma0=1}
if(missing(Num)){Num=1000000}
if(missing(a0)){a0=-1}
if(missing(b0)){b0=0.5}
if(missing(cutoff)){cutoff=0}
if(missing(lambda)){lambda=0}
if(missing(nu)){nu=Inf}
if(missing(regime)){stop('regimes need to be specified')}
if(length(regime)>2){stop('the length of regime should be either 1 or 2')}
if(missing(c_i)){b0=b0}
else
{
b0=uniroot(b0fun, c(-1, 1), tol = 0.0001, c_i=c_i, Sigma=Sigma, sigma1=sigma1, nu=nu, lambda=lambda, sigma0=sigma0)$root
}
if(missing(p_i)){a0=a0}
else
{
p_i=pifun(cutoff, a0, b0, Sigma, sigma0)
a0=uniroot(a0fun, c(-100, 100), tol = 0.0001, p_i=p_i, cutoff=cutoff, b0=b0, Sigma=Sigma, sigma0=sigma0)$root
}
m=dim(mu)[2]
z.b=qnorm(b);z.a.by.2=qnorm(1-a/2)
Sigma=CAR.Cov.Teeth(m, rho, tau)


#res:   matrix of dimension (#of treatment path X 1), where a row represents responding or non-responding that corresponds to a treatment path
#p_st2: matrix of dimension (#of treatment path X 1), the randomization probability matrix at stage-2, i.e. the rows represent the treatment path
#ga:    matrix of dimension (#of treatment path X 1), gamma, the response rates of initial treatments correspond to each row of res or p_st2


initr=as.matrix( rep(st1[,4], st1[,1]+st1[,2]) )
ga=as.matrix( rep(st1[,3], st1[,1]+st1[,2]) )
res=c();
p_st2=c()
for(i in 1:dim(st1)[1])
{
res=c( res, rep( c(1,0), st1[i,1:2]  ) )
p_st2=c( p_st2, rep( 1/st1[i,1:2], st1[i,1:2]  ) )
}
res=as.matrix(res);
p_st2=as.matrix(p_st2)

ga_comp=rep(0, length(regime))
mu_comp=array(0,c(2,m,length(regime)))
p2_comp=matrix(0, 2, length(regime))
p1_comp=rep(0,length(regime))

if(length(regime)==2)
{
ga_comp[1]=ga[ dtr[ regime[1],2 ], ]
ga_comp[2]=ga[ dtr[ regime[2],2 ], ]
mu_comp[,,1]=rbind( mu[ dtr[ regime[1],2 ], ], mu[ dtr[ regime[1],3 ], ] )
mu_comp[,,2]=rbind( mu[ dtr[ regime[2],2 ], ], mu[ dtr[ regime[2],3 ], ] )
p2_comp[,1]=rbind( p_st2[ dtr[ regime[1],2 ], ], p_st2[ dtr[ regime[1],3 ], ] )
p2_comp[,2]=rbind( p_st2[ dtr[ regime[2],2 ], ], p_st2[ dtr[ regime[2],3 ], ] )
p_st1=as.matrix( 
rep(
( colSums( t( 1/st1[,1:2] )*rbind( st1[,3],( 1-st1[,3] ) ) ) )^(-1) /sum( ( colSums( t( 1/st1[,1:2] )*rbind( st1[,3],( 1-st1[,3] ) ) ) )^(-1) ) 
, st1[,1]+st1[,2]
)
)
p1_comp[1]=p_st1[ dtr[ regime[1],2 ], ]
p1_comp[2]=p_st1[ dtr[ regime[2],2 ], ]
p1_comp=as.vector(p1_comp)

Yibard1R=MC.Var.Yibar.Mis(mu_comp[1,,1], Sigma, sigma1, lambda, nu, sigma0, Num, a0, b0, cutoff)
Sigd1R=Yibard1R$VarYi
mud1R=Yibard1R$mYi
Yibard1NR=MC.Var.Yibar.Mis(mu_comp[2,,1], Sigma, sigma1, lambda, nu, sigma0, Num, a0, b0, cutoff)
Sigd1NR=Yibard1NR$VarYi
mud1NR=Yibard1NR$mYi
Yibard2R=MC.Var.Yibar.Mis(mu_comp[1,,2], Sigma, sigma1, lambda, nu, sigma0, Num, a0, b0, cutoff)
Sigd2R=Yibard2R$VarYi
mud2R=Yibard2R$mYi
Yibard2NR=MC.Var.Yibar.Mis(mu_comp[2,,2], Sigma, sigma1, lambda, nu, sigma0, Num, a0, b0, cutoff)
Sigd2NR=Yibard2NR$VarYi
mud2NR=Yibard2NR$mYi
Del=abs( mud2R*ga_comp[2] + mud2NR*( 1-ga_comp[2] ) -( mud1R*ga_comp[1] + mud1NR*(1-ga_comp[1]) ) )
ybard2=mud2R*ga_comp[2] + mud2NR*( 1-ga_comp[2]) 
ybard1=mud1R*ga_comp[1] + mud1NR*(1-ga_comp[1])

vr1=( ga_comp[1]/( p1_comp[1]*p2_comp[1,1] ) )*( ( Sigd1R )+( 1-p1_comp[1]*p2_comp[1,1] )*( mud1R )^2 )
vnr1=( ( 1-ga_comp[1] )/( p1_comp[1]*p2_comp[2,1] ) )*( ( Sigd1NR )+( 1-p1_comp[1]*p2_comp[2,1] )*( mud1NR )^2 )
vrnrdiff1=ga_comp[1]*( 1-ga_comp[1] )*( mud1R - mud1NR )^2
vr2=( ga_comp[2]/( p1_comp[2]*p2_comp[1,2] ) )*( ( Sigd2R )+( 1-p1_comp[2]*p2_comp[1,2] )*( mud2R )^2 )
vnr2=( ( 1-ga_comp[2] )/( p1_comp[2]*p2_comp[2,2] ) )*( ( Sigd2NR )+( 1-p1_comp[2]*p2_comp[2,2] )*( mud2NR )^2 )
vrnrdiff2=ga_comp[2]*( 1-ga_comp[2] )*( mud2R - mud2NR )^2

if(dtr[ regime[1],4 ]==dtr[ regime[2],4 ])
{
covr=( ga_comp[1]/( p1_comp[1]*p2_comp[1,1] ) )*( ( Sigd1R ) + ( mud1R )^2 ) - ( ga_comp[1]*ga_comp[2] )*( mud1R * mud2R ) - ( ga_comp[1]*(1-ga_comp[2]) )*( mud1R * mud2NR ) - ( (1-ga_comp[1])*ga_comp[2] )*( mud1NR * mud2R ) - ( (1-ga_comp[1])*(1-ga_comp[2]) )*( mud1NR * mud2NR )                          
}
else if(dtr[ regime[1],4 ]!=dtr[ regime[2],4 ])
{
covr= - ( ga_comp[1]*ga_comp[2] )*( mud1R * mud2R ) - ( ga_comp[1]*(1-ga_comp[2]) )*( mud1R * mud2NR ) - ( (1-ga_comp[1])*ga_comp[2] )*( mud1NR * mud2R ) - ( (1-ga_comp[1])*(1-ga_comp[2]) )*( mud1NR * mud2NR )    
}
else{covr=0}
}
else if(length(regime)==1)
{
ga_comp[1]=ga[ dtr[ regime[1],2 ], ]
mu_comp[,,1]=rbind( mu[ dtr[ regime[1],2 ], ], mu[ dtr[ regime[1],3 ], ] )
p2_comp[,1]=rbind( p_st2[ dtr[ regime[1],2 ], ], p_st2[ dtr[ regime[1],3 ], ] )
p_st1=as.matrix( 
rep(
( colSums( t( 1/st1[,1:2] )*rbind( st1[,3],( 1-st1[,3] ) ) ) )^(-1) /sum( ( colSums( t( 1/st1[,1:2] )*rbind( st1[,3],( 1-st1[,3] ) ) ) )^(-1) ) 
, st1[,1]+st1[,2]
)
)
p1_comp[1]=p_st1[ dtr[ regime[1],2 ], ]
p1_comp=as.vector(p1_comp)


Yibard1R=MC.Var.Yibar.Mis(mu_comp[1,,1], Sigma, sigma1, lambda, nu, sigma0, Num, a0, b0, cutoff)
Sigd1R=Yibard1R$VarYi
mud1R=Yibard1R$mYi
Yibard1NR=MC.Var.Yibar.Mis(mu_comp[2,,1], Sigma, sigma1, lambda, nu, sigma0, Num, a0, b0, cutoff)
Sigd1NR=Yibard1NR$VarYi
mud1NR=Yibard1NR$mYi
Yibard2R=0
Sigd2R=0
mud2R=0
Yibard2NR=0
Sigd2NR=0
mud2NR=0
Del=abs( ( mud1R*ga_comp[1] + mud1NR*(1-ga_comp[1]) ) )
ybard2=0
ybard1=mud1R*ga_comp[1] + mud1NR*(1-ga_comp[1])

vr1=( ga_comp[1]/( p1_comp[1]*p2_comp[1,1] ) )*( ( Sigd1R )+( 1-p1_comp[1]*p2_comp[1,1] )*( mud1R )^2 )
vnr1=( ( 1-ga_comp[1] )/( p1_comp[1]*p2_comp[2,1] ) )*( ( Sigd1NR )+( 1-p1_comp[1]*p2_comp[2,1] )*( mud1NR )^2 )
vrnrdiff1=ga_comp[1]*( 1-ga_comp[1] )*( mud1R - mud1NR )^2
vr2=0
vnr2=0
vrnrdiff2=0
covr=0
}

sig.d1.sq=vr1+vnr1+vrnrdiff1;
sig.d2.sq=vr2+vnr2+vrnrdiff2;
sig.d1d2=covr;
sig.e.sq=sig.d1.sq+sig.d2.sq-2*sig.d1d2
N=( ( z.a.by.2-z.b )^2 )*( sig.e.sq/Del^2 )
Del_std=Del/sqrt(sig.e.sq/2)
return(list(
N=N,
sig.d1.sq=sig.d1.sq,
sig.d2.sq=sig.d2.sq,
sig.d1d2=sig.d1d2,
sig.e.sq=sig.e.sq,
Del=Del,
Del_std=Del_std,
ybard1=ybard1,
ybard2=ybard2,
initr=initr,
ga=ga,
res=res,
p_st2=p_st2,
p_st1=p_st1,
Sigma=Sigma
)
)
}
                                                 










                                             