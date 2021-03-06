% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ICV.R
\name{ICV}
\alias{ICV}
\title{The ICV function.}
\usage{
ICV(h, x)
}
\arguments{
\item{h}{numerical vector of bandwidth values (in the final scale),}

\item{x}{numerical vecror of data.}
}
\value{
The value of \eqn{ICV(h)} for given \eqn{h} and data (\eqn{x}).
}
\description{
Computing \eqn{ICV(h)}, the value of the ICV function, at a given bandwidth \eqn{h} (vector) for a data set \eqn{x} of size \eqn{n<12,058}. See Savchuk, Hart, and Sheather (2010).
}
\details{
Computation of  \eqn{ICV(h)} for given \eqn{h} (bandwidth vector) and \eqn{x} (data vector). The sample size \eqn{n<12,058}. The Gaussian kernel is to be used for computing the ultimate kernel density estimate. The parameters of the selection kernel are \eqn{(\alpha,\sigma)=(2.42, 5.06)}. The ICV bandwidth \code{\link{h_ICV}} is the minimizer of the ICV function.
}
\examples{
#Example 1. Computation of ICV(h) at h=0.4 for a random sample of size n=100 from a N(0,1) distribution.
ICV(0.4,rnorm(100))

#Example 2. (Calculations for a random sample of size n=250 from the separated bimodal density).
w=c(1/2,1/2)
mu=c(-3/2,3/2)
sdev=c(1/2,1/2)
# Generating a sample of size n=250 from the separated bimodal density of Marron and Wand (1992).
dat=mixnorm(250,w,mu,sdev)
h_OS=(243/(35*length(dat)*2*sqrt(pi)))^0.2*sd(dat)    # Computing the oversmoothed bandwidth.
h_opt=round(h_ICV(dat),digits=4)
harg=seq(0.1,3,len=100)
X11()
plot(harg,ICV(harg,x=dat),'l',lwd=3,xlab="h",ylab="ICV",cex.lab=1.7,cex.axis=1.7)
title(main="ICV(h)",cex.main=1.7)
lines(c(h_OS,h_OS),c(-0.5,0.5),lty="dashed",lwd=3)
legend(0.75,-0.05,legend="Vertical line shows the oversmothed bandwidth")
legend(1.35,0.1,legend=paste("h_ICV=",h_opt),cex=2,bty="n")
# Notice that the scale of the ICV function is such that its minimizer is the ICV bandwidth h_ICV.
# Thus, no additional rescaling of the ICV function's minimizer is needed to obtain the ICV bandwidth.
X11()
dens=density(dat,bw=h_opt)
plot(dens,main="",cex.lab=1.7,cex.axis=1.7,lwd=3,xlab=paste("n=250, ","h=h_ICV=",h_opt),ylim=c(0,0.45))
title(main="KDE based on h_ICV",cex.main=1.7)
arg=seq(-3.5,3.5,len=1000)
lines(arg,w[1]*dnorm(arg,mu[1],sd=sdev[1])+w[2]*dnorm(arg,mu[2],sd=sdev[2]),lwd=3,lty="dashed")
legend(-1,0.45,lty=c("solid","dashed"),lwd=c(3,3),legend=c("ICV estimate","True density"),bty="n")
}
\references{
Savchuk, O.Y., Hart, J.D., Sheather, S.J. (2010). Indirect cross-validation for density estimation. \emph{Journal of the American Statistical Association}, 105(489), 415-423.
}
\seealso{
\code{\link{h_ICV}}, \code{\link{C_ICV}}, \code{\link{L_ICV}}, \code{\link{MISE_mixnorm}}, \code{\link{KDE_ICV}}, \code{\link{LocICV}}.
}

