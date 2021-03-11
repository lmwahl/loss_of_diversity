## 
## Plots for environment change and extinction
##  for model with multiple haplotypes and environments
##   - haplotype frequency p vs. P_stand(p) 
##   - (cumulative) env change probability vs. haploid frequencies
##
## Make the plot by calling the following in R:
##
##  > source("env-change-extinct.R")
##  > plot.haplotypes.env() 



## Probability of evo rescue - via standing variation only 
## Pstand as a function of p
##  - default parameters given in arguments
Pstand <- function(p, N=10^3, sr=0.06, delta=0.01){
   1 - exp(-2*(sr-delta)*N*p)
}


## Plot quantities for multi-environment/haplotype model
plot.haplotypes.env <- function() {
   pdf("haplotypes-environments.pdf", height=2.5
     , title="Haplotypes and environments")
   par(mfrow=c(1,3), mar=c(4,4,2,1), cex=0.7)

   ## First plot: Pstand as a function of p 
   p <- seq(0,1, length.out=1000)
   plot(p, Pstand(p), t="l", xlab="Haplotype frequency, p"
      , ylab=expression(P[stand](p))
      , main =  "Probability of rescue", cex.main=0.8)
   lines(p,p, lty=3)
   mtext("(A)", 3, line=0.75, adj=-0.2, cex=0.75) 
   
   boxcol="grey80"  # colour of the rectangles 
   
   ## Second: "random" - no particular relationship between p and epsilon
   p <- c(0.2, 0.5, 0.2, 0.1); eps <- c(0.2, 0.3, 0.1, 0.4)
   cp <- cumsum(p)/sum(p)     # cumulative frequencies 
   ce <- cumsum(eps)/sum(eps) # cumulative env change prob 
   getprobs(p,eps)   # report the probs of rescue 
   
   plot(1, t="n", xlim=c(0,1), ylim=c(0,1)
      , xlab="Environment change probabilities", ylab="Haploid frequencies"
      , main="No particular relationship", cex.main=0.8)
   polygon(x=c(0,0, eps[1],eps[1]),
           y=c(0,p[1], p[1], 0), col=boxcol )
   polygon(x=c(eps[1],eps[1], ce[2],ce[2]),
           y=c(p[1], cp[2], cp[2], p[1]), col=boxcol )
   polygon(x=c(ce[2],ce[2], ce[3], ce[3]),
           y=c(cp[2], cp[3], cp[3], cp[2]), col=boxcol )
   polygon(x=c(ce[3],ce[3], ce[4], ce[4]),
           y=c(cp[3], cp[4], cp[4], cp[3]), col=boxcol )
   mtext("(B)", 3, line=0.75, adj=-0.2, cex=0.75) 
   
   ## Third: a mismatch between p and epsilon 
   p <- c(0.9, 0.05, 0.03, 0.02); eps <- c(0, 0.02, 0.04, 0.94)
   cp <- cumsum(p)/sum(p)
   ce <- cumsum(eps)/sum(eps)
   getprobs(p,eps)
   
   plot(1, t="n", xlim=c(0,1), ylim=c(0,1)
      , xlab="Environment change probabilities", ylab="Haploid frequencies"
      , main="Mismatch" , cex.main=0.8 )
   polygon(x=c(0,0, eps[1],eps[1]),
           y=c(0,p[1], p[1], 0), col=boxcol )
   polygon(x=c(eps[1],eps[1], ce[2],ce[2]),
           y=c(p[1], cp[2], cp[2], p[1]), col=boxcol )
   polygon(x=c(ce[2],ce[2], ce[3], ce[3]),
           y=c(cp[2], cp[3], cp[3], cp[2]), col=boxcol )
   polygon(x=c(ce[3],ce[3], ce[4], ce[4]),
           y=c(cp[3], cp[4], cp[4], cp[3]), col=boxcol )
   mtext("(C)", 3, line=0.75, adj=-0.2, cex=0.75) 
   
   dev.off()
}


## Compute and print probabilities of rescue after sweep/no-sweep
getprobs <- function(p,eps) {
   sweep.rescue <- sum(p*eps)
   no.sweep.rescue <- sum(Pstand(p)*eps)
   cat("P(rescue after no sweep)", no.sweep.rescue, "\n")
   cat("P(rescue after sweep)", sweep.rescue, "\n")
}



