
# a few first tests: look at the differences in the calculated discrimination values 
# using the simple (linear) and more comprehensive model

# variant 1 is supposed to represent deciduous broadleaf trees
# variant 2 is supposed to represent needle-leaf evergreen trees

Ca <- c(300:400)
Ci1 <- 0.7*Ca
Ci2 <- 0.6*Ca
An1 <- seq(7,1.2*7,length.out=101)
An2 <- seq(4,1.2*4,length.out=101)

#gm1 <- 0.175
#gm2 <- 0.078
gm1 <- 0.20
gm2 <- 0.10


Gamma_star <- 37.43     # Gamma* from Bernacchi et al. 2002 (note that values differ depending on whether gm is accounted for or not)
Cc1 <- Ci1 - An1/gm1
Cc2 <- Ci2 - An2/gm2

gs1 <- An1 / (Ca - Ci1)
gs2 <- An2 / (Ca - Ci2)


# values from Ubierna & Farquhar 2014 PCE, and in permil
a       <- 4.4
b       <- 30
b_prime <- 27
am      <- 1.8
f       <- 12


## linear model
delta_lin1 <- a + (b_prime - a) * Ci1/Ca
delta_lin2 <- a + (b_prime - a) * Ci2/Ca

## classical model (considering gm and photorespiration (Gamma*))
delta_clas1 <- a*(Ca-Ci1)/Ca + am*(Ci1-Cc1)/Ca + b*(Cc1/Ca) - f*(Gamma_star/Ca)
delta_clas2 <- a*(Ca-Ci2)/Ca + am*(Ci2-Cc2)/Ca + b*(Cc2/Ca) - f*(Gamma_star/Ca)






#############
### Plots ###
#############
## a) changes in discrimination
graphics.off()
lwd       <- 2
col_lin1  <- "darkgreen"
col_lin2  <- "brown"
col_clas1 <- "green2"
col_clas2 <- "orange"
par(mar=c(4,4.5,1.5,1.5))

plot(delta_lin1 ~ Ca,col=col_lin1,type="l",ylim=c(8.5,23),las=1,ylab=expression(Delta),tcl=-0.2,mgp=c(2.4,0.5,0),
     xlab=expression("C"[a]~"("*mu*"mol mol"^{-1}*")"),cex.lab=1.3,lwd=lwd)
points(delta_lin2 ~ Ca,col=col_lin2,type="l",lwd=lwd)
points(delta_clas1 ~ Ca,col=col_clas1,type="l",lwd=lwd)
points(delta_clas2 ~ Ca,col=col_clas2,type="l",lwd=lwd)
legend(x=330,y=13.6,legend=c("simple DBF","simple ENF","classical DBF","classical ENF"),
       col=c(col_lin1,col_lin2,col_clas1,col_clas2),lty=1,bty="n",
       x.intersp=0.6,seg.len=1,y.intersp=1.15,lwd=lwd)
mtext(side=3,at=305,line=-1.32,"(a)",cex=1.1)


dev.copy2pdf(file="C:/Profiles/jknauer/Desktop/13C_models_discrimination.pdf",width=7.3,height=4)




# calc_fract <- function(delta_atm,Discrimination,R_standard=0.0112372){
#   delta_atm      <- delta_atm/1000
#   Discrimination <- Discrimination/1000
#   
#   R_atm         <- (delta_atm + 1) * R_standard 
#   R_plant       <- R_atm / pmax(1e-12,(Discrimination + 1))
#   fractionation <- 1 / (1 + 1 / R_plant)
#   
#   return(fractionation)
# }
