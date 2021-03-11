## 
## Plot heterozygosity data from Spielman et al 
##
## Make the plot by calling the following in R:
## 
##  > source("het-plot.R") 
##  > plot.spielman.segments() 


## =========================================================
## Read in the data from a file 
d <- read.table("spielman-data-only-het-pairs.csv", sep=",", header=TRUE)
Hnl  <- na.omit(d$Hnl) # remove na (empty rows) 
Hl  <- na.omit(d$Hl)   # remove na (empty rows) 


## =========================================================
## Plot segments showing difference in H 
plot.spielman.segments <- function(){
   pdf("segments-spielman-data.pdf", height=2.3,width=6,
       title="Spielman H differences")
   par(mar=c(2,4,1,1),cex=0.65)
   ## Prepare index vector etc 
   ind <- 1:length(Hl)
   cola <- rep(2,length(Hl))  # first make everything red 
   cola[which(Hl>Hnl)]=1      # black if listed has higher H

   ## data sorted
   temp <- data.frame(Hnl,Hl,cola)
   ds <- temp[order(temp$Hnl),]

   ## Set up empty plot
   plot(ind,Hnl, t="n",
        ylab="Heterozygosity", xlab="",
        ylim=c(0,1)
      , bty="n", xaxt="n" )
   abline(h=0)
   mtext(expression(paste("Pairs of taxa (ordered by ", H[nl],")" ) )
       , side=1, line=1, cex=0.7) 

   ## Plot sorted segments
   segments(ind,ds$Hnl,ind,ds$Hl,  col=ds$cola)
   legend("topleft", c(expression(H[nl]>H[l]), expression(H[nl]<H[l])),
       lty=1, col=c(2,1), bty="n" )

   dev.off()
}


## =========================================================
## Plot density estimates and histogram of differences
plot.hist.density <- function(){ 
   pdf("hist-density-spielman-data.pdf",height=3,
       title="Heterozygosity distributions")
   op <- par(mfrow=c(1,2), mar=c(4,4,1,1), cex=0.7)

   ## First plot density estimates of H 
   dl <- density(Hl, from=0, to=1)
   peak.l <- dl$x[dl$y==max(dl$y)]
   dnl <- density(Hnl, from=0, to=1)
   peak.nl <- dnl$x[dnl$y==max(dnl$y)]

   plot(dl, main = "", xlab="Heterozygosity", col=2) 
   lines(dnl,  col=4)
   
   abline(v=c(peak.l, peak.nl), col=c(2,4), lty=2 )
   legend("topright",c("Non-threatened","Threatened"),
          lty=1, col=c(4,2))
   
   ## Next the histogram of differences: 
   hist(Hnl-Hl, breaks=25, main="", xlab="Differences: Hnl-Hl")
   dev.off()

   cat("The left peaks are at ", peak.l, " and ", peak.nl, "\n")
}
