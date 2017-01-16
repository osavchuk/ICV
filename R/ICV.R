#' The ICV bandwidth.
#'
#' Calculation of the ICV bandwidth for the Gaussian density estimator corresponding to expression (12) of Savchuk, Hart, and Sheather (2010).
#'
#' Computing the ICV bandwidth for a univariate numerical data set of size \eqn{n<12,058}. The ICV bandwidth is consistent for the MISE optimal bandwidth (see Wand and Jones (1995)). The Gaussian kernel is used for computing the ultimate density estimate. The following values of the paramaters of the selection kernel \code{\link{L_ICV}} are used: \eqn{(\alpha,\sigma)=(2.42, 5.06)}. The ICV bandwidth does not exceed the oversmoothed bandwidth of Terrell (1990). See expression (12) of Savchuk et al. (2010).
#' @param x numerical vector of data.

#' @return The ICV bandwidth.
#' @references
#' \itemize{
#'   \item Savchuk, O.Y., Hart, J.D., Sheather, S.J. (2010). Indirect cross-validation for density estimation. \emph{Journal of the American Statistical Association}, 105(489), 415-423.
#'   \item Wand, M.P. and Jones, M.C. (1995). Kernel Smoothing. Chapman and Hall, London.
#'   \item Terrel, G. (1990). The maximum smoothing principle in density estimation. \emph{Journal of the American Statistical Association}, 85, 470-477.
#' }
#' @seealso \code{\link{ICV}}, \code{\link{C_ICV}}, \code{\link{L_ICV}}, \code{\link{MISE_mixnorm}}, \code{\link{KDE_ICV}}, \code{\link{LocICV}}.
#' @examples
#' # ICV bandwidth for a random sample of size n=100 from a N(0,1) density.
#' h_ICV(rnorm(100))
h_ICV=function(x){
  # x is the data vector
  h_OS=(243/(35*length(x)*2*sqrt(pi)))^0.2*sd(x)      # computing the oversmoothed bandwidth
  optimize(ICV,c(0.001,h_OS),tol=0.001,x=x)$minimum
}


#' The ICV function.
#'
#' Computing \eqn{ICV(h)}, the value of the ICV function, at a given bandwidth \eqn{h} (vector) for a data set \eqn{x} of size \eqn{n<12,058}. See Savchuk, Hart, and Sheather (2010).
#'
#' Computation of  \eqn{ICV(h)} for given \eqn{h} (bandwidth vector) and \eqn{x} (data vector). The sample size \eqn{n<12,058}. The Gaussian kernel is to be used for computing the ultimate kernel density estimate. The parameters of the selection kernel are \eqn{(\alpha,\sigma)=(2.42, 5.06)}. The ICV bandwidth \code{\link{h_ICV}} is the minimizer of the ICV function.
#' @param h numerical vector of bandwidth values (in the final scale),
#' @param x numerical vecror of data.
#' @return The value of \eqn{ICV(h)} for given \eqn{h} and data (\eqn{x}).
#' @references Savchuk, O.Y., Hart, J.D., Sheather, S.J. (2010). Indirect cross-validation for density estimation. \emph{Journal of the American Statistical Association}, 105(489), 415-423.
#' @seealso \code{\link{h_ICV}}, \code{\link{C_ICV}}, \code{\link{L_ICV}}, \code{\link{MISE_mixnorm}}, \code{\link{KDE_ICV}}, \code{\link{LocICV}}.
#' @examples
#' #Example 1. Computation of ICV(h) at h=0.4 for a random sample of size n=100 from a N(0,1) distribution.
#' ICV(0.4,rnorm(100))
#'
#' #Example 2. (Calculations for a random sample of size n=250 from the separated bimodal density).
#' w=c(1/2,1/2)
#' mu=c(-3/2,3/2)
#' sdev=c(1/2,1/2)
#' # Generating a sample of size n=250 from the separated bimodal density of Marron and Wand (1992).
#' dat=mixnorm(250,w,mu,sdev)
#' h_OS=(243/(35*length(dat)*2*sqrt(pi)))^0.2*sd(dat)    # Computing the oversmoothed bandwidth.
#' h_opt=round(h_ICV(dat),digits=4)
#' harg=seq(0.1,3,len=100)
#' X11()
#' plot(harg,ICV(harg,x=dat),'l',lwd=3,xlab="h",ylab="ICV",cex.lab=1.7,cex.axis=1.7)
#' title(main="ICV(h)",cex.main=1.7)
#' lines(c(h_OS,h_OS),c(-0.5,0.5),lty="dashed",lwd=3)
#' legend(0.75,-0.05,legend="Vertical line shows the oversmothed bandwidth")
#' legend(1.35,0.1,legend=paste("h_ICV=",h_opt),cex=2,bty="n")
#' # Notice that the scale of the ICV function is such that its minimizer is the ICV bandwidth h_ICV.
#' # Thus, no additional rescaling of the ICV function's minimizer is needed to obtain the ICV bandwidth.
#' X11()
#' dens=density(dat,bw=h_opt)
#' plot(dens,main="",cex.lab=1.7,cex.axis=1.7,lwd=3,xlab=paste("n=250, ","h=h_ICV=",h_opt),ylim=c(0,0.45))
#' title(main="KDE based on h_ICV",cex.main=1.7)
#' arg=seq(-3.5,3.5,len=1000)
#' lines(arg,w[1]*dnorm(arg,mu[1],sd=sdev[1])+w[2]*dnorm(arg,mu[2],sd=sdev[2]),lwd=3,lty="dashed")
#' legend(-1,0.45,lty=c("solid","dashed"),lwd=c(3,3),legend=c("ICV estimate","True density"),bty="n")
ICV=function(h,x){       #ICV function from JASA paper
    # h is a bandwidth(vector) in the final scale;
    # x is a data vecror
    # The choice of (alpha,sigma)=(2.42,5.06) is used for n<12,058.
    ngrid=length(h)
    a=2.4233
    s=5.06
    C_ICV=3.33496
    h=h/C_ICV
    n=length(x)
    ICV=1:ngrid
    arg=matrix(x,n,1)%*%matrix(1,1,n)-matrix(1,n,1)%*%matrix(x,1,n)
    for(i in 1:ngrid){
      S1=((1+a)^2*dnorm(arg,sd=h[i]*sqrt(2))-2*a*(1+a)*dnorm(arg,sd=h[i]*sqrt(1+s^2))+a^2*dnorm(arg,sd=h[i]*s*sqrt(2)))/n^2
      S2=(2*a*dnorm(arg,sd=h[i]*s)-2*(1+a)*dnorm(arg,sd=h[i]))/(n*(n-1))
      diag=2*a/((n-1)*h[i]*s*sqrt(2*pi))-2*(1+a)/((n-1)*h[i]*sqrt(2*pi))
      ICV[i]=sum(S1+S2)-diag
    }
    h=h*C_ICV
    ICV
 }



#' Generating a random sample from the specified mixture of normal distributions.
#'
#' Generating a random sample of size \eqn{n} from the normal mixture defined by  expression (2.3) of Marron and Wand (1992).
#'
#' Producing a random sample of size \eqn{n} from the normal mixture defined by the vector of weights \eqn{w}, the vector of means \eqn{\mu}, and the vector of standard deviations \eqn{\sigma}. See Marron and Wand (1992). It is assumed that the normals are defined as parsimonious as possible. The normal distributions in the mixture should be ordered such that the means in \eqn{\mu} are arranged in a nondecreasing order.
#' @param  n  desired sample size,
#' @param  w vector of weighs (positive numbers between 0 and 1 that add up to one),
#' @param  mu vector of means,
#' @param  sdev vector of standard deviations.
#' @return A random sample of size \eqn{n} from the specified mixture of normals.
#' @references Marron, J.S., Wand, M.P. (1992). Exact Mean Integrated Squared Error. \emph{The Annals of Statistics}, 20(2), 712-736.
#' @seealso \code{\link{ISE_mixnorm}}, \code{\link{h_isemixnorm}}, \code{\link{MISE_mixnorm}}.
#' @examples
#' # Generating a sample of size n=300 from the separated bimodal density of Marron and Wand (1992).
#' w=c(0.5,0.5)
#' mu=c(-3/2,3/2)
#' sdev=c(1/2,1/2)
#' dat=mixnorm(300,w,mu,sdev) # generated data vector
#' arg=seq(-4,4,len=1000)  # argument
#' f=w[1]*dnorm(arg,mu[1],sd=sdev[1])+w[2]*dnorm(arg,mu[2],sd=sdev[2])     # true density
#' X11()
#' hist(dat,freq=F,ylab="",main="",cex.lab=1.7,cex.axis=1.7,xlim=c(-4,4),lwd=2,ylim=c(0,0.45),col='grey')
#' title(main="Separated bimodal density",cex.main=1.7)
#' legend(-5,0.4,legend="n=300",cex=2,bty="n")
#' lines(arg,f,lwd=3,'l')
mixnorm=function(n,w,mu,sdev){
  if((sum(w)<0.99)|(sum(w)>1.01)) print("The weights in w do not add up to 1")
  else{
    x=rnorm(n,mean=mu[1],sd=sdev[1])
    if(length(w)!=1){
      u=runif(n)
      for(i in 1:length(w)){
        f=(1:n)[(u>=sum(w[1:(i-1)]))&(u<sum(w[1:i]))]
        x[f]=rnorm(length(f),mean=mu[i],sd=sdev[i])
      }
    }
    x
  }
}


#'  The ISE function in the case when the underlying density is the specified mixture of normal distributions.
#'
#' Computing \eqn{ISE(h)} for given \eqn{h} in the case when the underlying density is the specified mixture of normal distributions and the Gaussian kernel is used to compute the ultimate density estimate.
#'
#' Computing \eqn{ISE(h)} in the case when the true density is the mixture of normal distributions defined by the vector of weights \eqn{w}, the vector of means \eqn{\mu}, and the vector of standard deviations \eqn{\sigma}. See expression (2.3) of Marron and Wand (1992). It is assumed that the normals are defined as parsimonious as possible. The normal distributions in the mixture should be ordered such that the means in \eqn{\mu} are arranged in a nondecreasing order. The Gaussian kernel is to be used for computing the ultimate density estimate.
#' @param  h  numerical vector of bandwidth values,
#' @param  x numerical vector of data,
#' @param  w vector of weighs (positive numbers between 0 and 1 that add up to one),
#' @param  mu vector of means,
#' @param  sdev vector of standard deviations.
#' @return The vector of ISE values corresponding to the specifies values of \eqn{h}.
#' @references Marron, J.S., Wand, M.P. (1992). Exact Mean Integrated Squared Error. \emph{The Annals of Statistics}, 20(2), 712-736.
#' @seealso \code{\link{mixnorm}}, \code{\link{h_isemixnorm}}, \code{\link{MISE_mixnorm}}.
#' @examples
#' harg=seq(0.01,1,len=100)
#' w=c(3/4,1/4)
#' mu=c(0,3/2)
#' sdev=c(1,1/3)
#' # The vectors w, mu, and sdev define the skewed bimodal density of Marron and Wand (1992).
#' dat=mixnorm(300,w,mu,sdev)
#' h_ISE=round(h_isemixnorm(dat,w,mu,sdev),digits=4)
#' ISEarray=ISE_mixnorm(harg,dat,w,mu,sdev)
#' X11()
#' plot(harg,ISEarray,'l',lwd=3,xlab="h, n=300",ylab="ISE",cex.lab=1.7,cex.axis=1.7,main="")
#' title(main="ISE(h)",cex.main=1.7)
#' legend(0.2,0.08,legend=paste("h_ISE=",h_ISE),cex=2)
ISE_mixnorm=function(h,x,w,mu,sdev){
   n=length(x)
   k=length(w)
   if((sum(w)<0.99)|(sum(w)>1.01)) print("The weights in w do not add up to 1")
   else{
     ngrid=length(h)
     arg=matrix(x,n,1)%*%matrix(1,1,n)-matrix(1,n,1)%*%matrix(x,1,n)
     m=matrix(mu,k,1)%*%matrix(1,1,k)-matrix(1,k,1)%*%matrix(mu,1,k)
     var=matrix(sdev^2,k,1)%*%matrix(1,1,k)+matrix(1,k,1)%*%matrix(sdev^2,1,k)
     fi=dnorm(m,sd=sqrt(var))
     R.f=sum(diag(w)%*%fi%*%diag(w))
     dif=matrix(x,n,1)%*%matrix(1,1,k)-matrix(1,n,1)%*%matrix(mu,1,k)
     ISE=1:ngrid
     for(i in 1:ngrid){
       R.fhat=1/n^2*sum(dnorm(arg,sd=h[i]*sqrt(2)))
       SD=matrix(1,n,1)%*%matrix(sqrt(sdev^2+h[i]^2),1,k)
       cross=sum(dnorm(dif,sd=SD)%*%matrix(w,k,1))/n
       ISE[i]=R.fhat+R.f-2*cross
     }
     ISE
  }
}

#' The ISE-optimal bandwidth in the case when the true density is the specified mixture of normal distributions.
#'
#' Computing the ISE-optimal bandwidth in the case when the true density is the specified mixture of normal distributions and the Gaussian kernel is used to compute the ultimate density estimate.
#'
#' Computing the ISE-optimal bandwidth (i.e. the minimizer of the ISE function) in the case when the true density is the mixture of normal distributions defined by the vector of weights \eqn{w}, the vector of means \eqn{\mu}, and the vector of standard deviations \eqn{\sigma}. See expression (2.3) of Marron and Wand (1992). It is assumed that the normals are defined as parsimonious as possible. The normal distributions in the mixture should be ordered such that the means in \eqn{\mu} are sorted in a nondecreasing order. The Gaussian kernel is used for computing the ultimate density estimate.
#' @param  x numerical vector of data,
#' @param  w vector of weighs (positive numbers between 0 and 1 that add up to one),
#' @param  mu vector of means,
#' @param  sdev vector of standard deviations.
#' @return The ISE-optimal bandwidth.
#' @references Marron, J.S., Wand, M.P. (1992). Exact Mean Integrated Squared Error. \emph{The Annals of Statistics}, 20(2), 712-736.
#' @seealso \code{\link{mixnorm}}, \code{\link{ISE_mixnorm}}, \code{\link{MISE_mixnorm}}.
#' @examples
#' # ISE optimal bandwidth for a random sample of size n=100 generated from a normal mixture defined by
#' # w=c(1/5,1/5,3/5), mu=(0,1/2,13/12), sdev=c(1,2/3,5/9).
#' # This corresponds to the skewed unimodal density of Marron and Wand (1992).
#' h_isemixnorm(rnorm(100),c(1/5,1/5,3/5),c(0,1/2,13/12),c(1,2/3,5/9))
h_isemixnorm=function(x,w,mu,sdev){
  if((sum(w)<0.99)|(sum(w)>1.01)) print("The weights in w do not add up to 1")
  else
    optimize(ISE_mixnorm,c(0.01,2),tol=0.001,x=x,w=w,mu=mu,sdev=sdev)$minimum
}



#'  The MISE function in the case when the true density is the specified mixture of normal distributions and the selection kernel  \code{\link{L_ICV}} is used in the cross-validation stage.
#'
#' Computing \eqn{MISE(h)} for given \eqn{h} in the case when the true density is the specified mixture of normal distributions and the kernel is \code{\link{L_ICV}} defined by expression (4) of Savchuk, Hart, and Sheather (2010).
#'
#' Calculation of \eqn{MISE(h)} in the case when the true density is the mixture of normal distributions defined by the vector of weights \eqn{w}, the vector of means \eqn{\mu}, and the vector of standard deviations \eqn{\sigma}. See expression (2.3) of Marron and Wand (1992). It is assumed that the normals are defined as parsimonious as possible. The normal distributions in the mixture should be ordered such that the means in \eqn{\mu} are arranged in a nondecreasing order. The MISE function is based on the selection kernel \code{\link{L_ICV}} defined by expression (4) of Savchuk, Hart, and Sheather (2010). Notice that the Gaussian kernel \eqn{\phi} is the special case of  \code{\link{L_ICV}} given that \strong{(Case 1)} \eqn{\alpha=0}, \eqn{\sigma>0} or \strong{(Case 2)} \eqn{\sigma=1}, \eqn{-\infty<\alpha<\infty}.
#' @param  h  numerical vector of bandwidth values,
#' @param  n  sample size,
#' @param  alpha  first parameter of the selection kernel,
#' @param  sigma  second parameter of the selection kernel,
#' @param  w  vector of weighs (positive numbers between 0 and 1 that add up to one),
#' @param  mu vector of means,
#' @param  sdev vector of standard deviations.
#' @return The vector of MISE values corresponding to the specified values of \eqn{h}.
#' @references
#' \itemize{
#'   \item Savchuk, O.Y., Hart, J.D., Sheather, S.J. (2010). Indirect cross-validation for density estimation. \emph{Journal of the American Statistical Association}, 105(489), 415-423.
#'   \item  Marron, J.S., Wand, M.P. (1992). Exact Mean Integrated Squared Error. \emph{The Annals of Statistics}, 20(2), 712-736.
#' }
#' @seealso \code{\link{mixnorm}}, \code{\link{ISE_mixnorm}}, \code{\link{h_isemixnorm}}, \code{\link{L_ICV}}, \code{\link{ICV}}, \code{\link{h_ICV}},  \code{\link{C_ICV}}.
#' @examples
#' # Example 1. MISE for the separated bimodal density of Marron and Wand (1992).
#' # in the case (alpha,sigma)=(2.42,5.06), n=100.
#' harray=seq(0.05,1,len=1000)
#' w=c(1/2,1/2)
#' m=c(-3/2,3/2)
#' s=c(1/2,1/2)
#' MISEarray=MISE_mixnorm(harray,100,2.42,5.06,w,m,s)
#' hopt=round(optimize(MISE_mixnorm,c(0.01,1),n=100,alpha=2.42,sigma=5.06,w=w,mu=m,sdev=s)$minimum,digits=4)
#' X11()
#' plot(harray,MISEarray,'l',lwd=3,xlab="h",ylab="MISE",cex.lab=1.7,cex.axis=1.7,main="")
#' title(main="MISE(h) for the separated bimodal density",cex.main=1.5)
#' legend(0.45,0.45,legend=c(paste("h_MISE=",hopt),"n=100"),bty="n",cex=1.7)
#'
#' # Example 2. MISE for the N(0,1) density in the case of the Gaussian kernel and n=500.
#' harray=seq(0.03,1,len=1000)
#' MISEarray=MISE_mixnorm(harray,500,1,1,1,0,1)
#' hopt=round(optimize(MISE_mixnorm,c(0.01,1),n=500,alpha=1,sigma=1,w=1,mu=0,sdev=1)$minimum,digits=4)
#' X11()
#' plot(harray,MISEarray,'l',lwd=3,xlab="h",ylab="MISE",cex.lab=1.7,cex.axis=1.7,main="")
#' title(main="MISE(h) for the standard normal density",cex.main=1.7)
#' legend(0.2,0.02,legend=c(paste("h_MISE=",hopt),"n=500"),bty="n",cex=1.7)
MISE_mixnorm=function(h,n,alpha,sigma,w,mu,sdev){
  # h is a bandwidth (vector); n=sample size;
  # alpha and sigma are parameters of the selection kernel;
  # vectors w,mu and sigma define the mixture of normals density
  ngrid=length(h)
  if((sum(w)<0.99)|(sum(w)>1.01)) print("The weights in w do not add up to 1")
  else{
    MISE=1:ngrid
    k=length(w)
    arg=matrix(mu,k,1)%*%matrix(1,1,k)-matrix(1,k,1)%*%matrix(mu,1,k)
    var=matrix(sdev^2,k,1)%*%matrix(1,1,k)+matrix(1,k,1)%*%matrix(sdev^2,1,k)
    ct=(1+alpha)^2/(2*sqrt(pi))+alpha^2/(2*sqrt(pi)*sigma)-sqrt(2/pi)*alpha*(1+alpha)/sqrt(1+sigma^2)
    Rf=dnorm(arg,sd=sqrt(var))
    for(i in 1:ngrid){
      R=1/(n*h[i])*ct
      var1=var+2*h[i]^2
      var2=var+h[i]^2*(1+sigma^2)
      var3=var+2*h[i]^2*sigma^2
      var4=var+h[i]^2
      var5=var+h[i]^2*sigma^2
      bracket=(n-1)*(1+alpha)^2/n*dnorm(arg,sd=sqrt(var1))-2*(n-1)*alpha*(1+alpha)/n*dnorm(arg,sd=sqrt(var2))+(n-1)*alpha^2/n*dnorm(arg,sd=sqrt(var3))-2*(1+alpha)*dnorm(arg,sd=sqrt(var4))+2*alpha*dnorm(arg,sd=sqrt(var5))+Rf
      doublesum=sum(diag(w)%*%bracket%*%diag(w))
      MISE[i]=R+doublesum
    }
    MISE
  }
}
####################################################################################
#' The ICV rescaling constant.
#'
#' Computing the ICV rescaling constant defined by expression (3) of Savchuk, Hart, and Sheather (2010).
#'
#' Calculation of the ICV rescaling constant \eqn{C} defined by (3) in Savchuk, Hart, and Sheather (2010). The constant is a function of the parameters \eqn{(\alpha,\sigma)} of the selection kernel \code{\link{L_ICV}} defined by expression (4) in the same article. The Gaussian kernel is to be used for computing the ultimate density estimate.
#' @param alpha first parameter of the selection kernel,
#' @param sigma second parameter of the selection kernel.
#' @return The ICV rescaling constant \eqn{C}.
#' @references  Savchuk, O.Y., Hart, J.D., Sheather, S.J. (2010). Indirect cross-validation for density estimation. \emph{Journal of the American Statistical Association}, 105(489), 415-423.
#' @seealso \code{\link{ICV}}, \code{\link{h_ICV}}, \code{\link{L_ICV}}, \code{\link{MISE_mixnorm}}, \code{\link{KDE_ICV}}, \code{\link{LocICV}}.
#' @examples
#' # ICV rescaling constant for the selection kernel with (alpha,sigma)=(2.42,5.06).
#' C_ICV(2.42,5.06)
C_ICV=function(alpha,sigma){
  R.phi=1/(2*sqrt(pi))
  R.L=(1+alpha)^2/(2*sqrt(pi))+alpha^2/(2*sqrt(pi)*sigma)-sqrt(2/pi)*alpha*(1+alpha)/sqrt(1+sigma^2)
  m2L=(1+alpha-alpha*sigma^2)
 (R.phi*m2L^2/R.L)^0.2
}

#' The ICV selection kernel.
#'
#' The ICV selection kernel \eqn{L} defined by expression (4) of Savchuk, Hart, and Sheather (2010).
#'
#' The ICV selection kernel \eqn{L(u;\alpha,\sigma)=(1+\alpha)\phi(u)-\alpha\phi(u/\sigma)/\sigma}, where \eqn{\phi} denotes the Gaussian kernel.
#' @param u numerical argument of the selection kernel,
#' @param alpha first parameter of the selection kernel,
#' @param sigma second parameter of the selection kernel.
#' @return The value of \eqn{L(u;\alpha,\sigma)}.
#' @references  Savchuk, O.Y., Hart, J.D., Sheather, S.J. (2010). Indirect cross-validation for density estimation. \emph{Journal of the American Statistical Association}, 105(489), 415-423.
#' @seealso \code{\link{ICV}}, \code{\link{h_ICV}}, \code{\link{C_ICV}}, \code{\link{MISE_mixnorm}}, \code{\link{KDE_ICV}}, \code{\link{LocICV}}.
#' @examples
#' # Graph of the ICV selection kernel with (alpha,sigma)=(2.42,5.06).
#' u=seq(-10,10,len=1000)
#' kern=L_ICV(u,2.42,5.06)
#' X11()
#' plot(u,kern,'l',lwd=2,ylim=c(-0.2,1.2),ylab="kernel",cex.lab=1.7,cex.axis=1.7,main="")
#' lines(u,dnorm(u),lwd=3,lty="dashed")
#' title(main="Selection kernel with (alpha,sigma)=(2.42,5.06)",cex.main=1.6)
#' legend(-11, 1.2, legend=c("ICV kernel","Gaussian kernel"),lwd=c(3,3),lty=c(1,2),bty="n",cex=1.3)
L_ICV=function(u,alpha,sigma)
  (1+alpha)*dnorm(u)-alpha/sigma*dnorm(u/sigma)
  
  
#' The local ICV function.
#'
#' Computing the local ICV function at the given estimation point, as explained in Section 6 of Savchuk, Hart, and Sheather (2010).
#'
#' Calculation of the local ICV function at the given estimation point xest. The Gaussian kernel is used for local weighting. The ultimate kernel density estimate is computed based on the Gaussian kernel. The parameters of the selection kernel \code{\link{L_ICV}} are \eqn{\alpha} and \eqn{\sigma}. The minimizer of the local ICV function is to be used in computing the ultimate density estimate without additional rescaling. Parameter \eqn{\eta} is a smoothing parameter that determines the degree to which the cross-validation is local. A suggested value of \eqn{\eta} is \eqn{\eta=R/20}, where \eqn{R} is the range of data.
#' @param h bandwidth (scalar) in the final scale,
#' @param xest estimation point (scalar),
#' @param x numerical vector of data,
#' @param eta smoothing parameter,
#' @param alpha first parameter of the selection kernel,
#' @param sigma second parameter of the selection kernel.
#' @return The value of the local ICV function at the fixed estimation point and for the specified value of the bandwidth (see Section 6 of Savchuk, Hart, and Sheather (2010)).
#' @references
#' \itemize{
#'   \item Savchuk, O.Y., Hart, J.D., Sheather, S.J. (2010). Indirect cross-validation for density estimation. \emph{Journal of the American Statistical Association}, 105(489), 415-423.
#'   \item Savchuk, O.Y., Hart, J.D., Sheather, S.J. (2009). An empirical study of indirect cross-validation. \emph{Nonparametric Statistics and Mixture Models: A Festschrift in Honor of Thomas P. Hettmansperger.} World Scientific Publishing, 288-308.
#'   \item Hall, P., and Schukany, W. R. (1989). A local cross-validation algorithm. \emph{Statistics and Probability Letters}, 8, 109-117.
#' }
#' @seealso \code{\link{h_ICV}}, \code{\link{C_ICV}}, \code{\link{L_ICV}}, \code{\link{MISE_mixnorm}}, \code{\link{ICV}}, \code{\link{KDE_ICV}}.
#' @examples
#' # Local ICV function for a random sample of size n=150 from the kurtotic density of Marron and Wand (1992).
#' dat=mixnorm(150,c(2/3,1/3),c(0,0),c(1,1/10))
#' a=2.42   # alpha
#' s=5.06   # sigma
#' harg=seq(0.025,1,len=100)
#' Xest=0.1    # estimation point
#' LocICV_Xest=numeric(length(harg))
#' for(i in 1:length(harg))
#'   LocICV_Xest[i]=LocICV(harg[i],Xest,dat,0.2,a,s)
#' h_Xest=optimize(LocICV,c(0.001,0.2),tol=0.001,xest=Xest,eta=0.2,x=dat,alpha=a,sigma=s)$minimum
#' h_Xest=round(h_Xest,digits=4)
#' X11()
#' plot(harg,LocICV_Xest,'l',lwd=3,xlab="harg",ylab="LocICV_Xest",main="",cex.lab=1.7, cex.axis=1.7)
#' title(main=paste("Local ICV function at x=",Xest),cex.main=1.7)
#' legend(0.1,max(LocICV_Xest),legend=paste("h_x=",h_Xest),cex=1.7)
#' legend(0.2,max(LocICV_Xest)-0.15,legend="Note: first local minimizer is used", cex=1.5,bty="n")
LocICV=function(h,xest,x,eta,alpha,sigma){
  h=h/C_ICV(alpha,sigma)
  n=length(x)
  arg1=matrix(x,n,1)%*%matrix(1,1,n)-matrix(1,n,1)%*%matrix(x,1,n)
  arg2=matrix(x,n,1)%*%matrix(1,1,n)+matrix(1,n,1)%*%matrix(x,1,n)
  arg3=matrix(x,n,1)%*%matrix(1,1,n)+sigma^2*matrix(1,n,1)%*%matrix(x,1,n)
  M1=(1+alpha)^2/n^2*dnorm(arg1,sd=h*sqrt(2))*dnorm(arg2/2-xest,sd=sqrt(eta^2+h^2/2))
  M2=2*alpha*(1+alpha)/n^2*dnorm(arg1,sd=h*sqrt(1+sigma^2))*dnorm(arg3/(1+sigma^2)-xest,sd=sqrt(eta^2+h^2*sigma^2/(1+sigma^2)))
  M3=alpha^2/n^2*dnorm(arg1,sd=h*sigma*sqrt(2))*dnorm(arg2/2-xest,sd=sqrt(eta^2+h^2*sigma^2/2))
  M4=2*(1+alpha)/(n*(n-1))*dnorm(arg1,sd=h)*(matrix(dnorm(x-xest,sd=eta))%*%matrix(1,1,n))
  M5=2*alpha/(n*(n-1))*dnorm(arg1,sd=h*sigma)*(matrix(dnorm(x-xest,sd=eta))%*%matrix(1,1,n))
  M6=2*(1+alpha-alpha/sigma)/(n*(n-1)*h*sqrt(2*pi))*sum(dnorm(x-xest,sd=eta))
  sum(M1-M2+M3-M4+M5)+M6
}


#' Computing the kernel density estimate based on the ICV bandwidth.
#'
#' Computing the Gaussian density estimate based on \code{\link{h_ICV}}.
#'
#' Computing the Gaussian density estimate based on \code{\link{h_ICV}}. The ICV selection kernel \code{\link{L_ICV}} is based on \eqn{(\alpha,\sigma)=(2.42,5.06)}.
#' @param x numerical vector of data.
#' @return A list containing the following components:
#' \describe{
#' \item{\eqn{arg}}{vector of sorted values of the argument at which the density estmate is computed,}
#' \item{\eqn{y}}{vector of density estimates at the corresponding values of the argument.}
#' }
#' The function also produces a graph of the resulting ICV kernel density estimate.
#' @references
#' \itemize{
#'   \item Savchuk, O.Y., Hart, J.D., Sheather, S.J. (2010). Indirect cross-validation for density estimation. \emph{Journal of the American Statistical Association}, 105(489), 415-423.
#'   \item Savchuk, O.Y., Hart, J.D., Sheather, S.J. (2009). An empirical study of indirect cross-validation. \emph{Nonparametric Statistics and Mixture Models: A Festschrift in Honor of Thomas P. Hettmansperger.} World Scientific Publishing, 288-308.
#' }
#' @seealso \code{\link{ICV}}, \code{\link{h_ICV}}, \code{\link{L_ICV}}, \code{\link{LocICV}}, \code{\link{C_ICV}}.
#' @examples
#' #Example (Density estimate for eruption duration of the Old Faithful Geyser).
#' data=faithful[[1]]
#' dens=KDE_ICV(data)
KDE_ICV=function(x){
  h.icv=round(h_ICV(x),digits=4)
  d=density(x,bw=h.icv)
  arg=d[[1]]
  y=d[[2]]
  X11()
  plot(d,cex.lab=1.7,cex.axis=1.7,xlab=paste(" h_ICV=",h.icv),ylab="density",main="",lwd=3)
  title(main="KDE based on h_ICV",cex.main=1.7)
  list(arg,y)
}





