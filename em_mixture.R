em_mixture <- function(x, p, mu, sd, eps=1e-8, maxit=1000) {
  mu1 <- mu[1]
  mu2 <- mu[2]
  sigsqrd1 <- sd[1]
  sigsqrd2 <- sd[2]
  mx <- mean(x)
  const <- length(x) * 0.918938533204673 # i.e., times log(2*pi)/2
  dl <- 1 + eps
  iter<-0
  ll <- rep(0, maxit+1)
  a1<-(x-mu1)^2
  b1<-(p/sqrt(sigsqrd1))*exp(-a1/2/sigsqrd1)
  a2<-(x-mu2)^2
  b2<-((1-p)/sqrt(sigsqrd2))*exp(-a2/2/sigsqrd2)
  l <- sum(log(b1+b2))

  while (dl>eps && iter<maxit && !is.na(dl)) {
    iter<-iter+1
    ll[iter] <- l
    postprobs <- b1/(b1+b2)
    p<-mean(postprobs)
    mu1<-mean(postprobs*x)/p
    mu2<-(mx-p*mu1)/(1-p)
    if (arbvar)  {
      sigsqrd1<-mean(postprobs*a1)/p
      sigsqrd2<-mean((1-postprobs)*a2)/(1-p)
    } else {
      sigsqrd1 <- sigsqrd2 <- mean(postprobs*a1 + (1-postprobs)*a2)
    }
    a1<-(x-mu1)^2; b1<-(p/sqrt(sigsqrd1))*exp(-a1/2/sigsqrd1)
    a2<-(x-mu2)^2; b2<-((1-p)/sqrt(sigsqrd2))*exp(-a2/2/sigsqrd2)

    oldl<-l
    l <- sum(log(b1+b2))
    dl<-l-oldl
  }
  iter <- iter+1
  ll[iter] <- l
  postprobs <- cbind(postprobs, 1-postprobs)
  colnames(postprobs) <- c(paste("comp", ".", 1:2, sep = ""))
  out <- list(x=x, p = c(p,1-p), mu = c(mu1, mu2),
       sigma = sqrt(c(sigsqrd1, sigsqrd2)[1:(1+arbvar)]),
       loglik = l - const, posterior = postprobs,
       all.loglik=ll[1:iter] - const,
       restarts=0, ft="normalmixEM")
  class(out) <- "mixEM"
  out
}

