madata=read.csv("https://figshare.com/ndownloader/files/14460386")

library(metafor)
d2=escalc(data=madata, n1i = Ne, n2i = Nc, m1i = Me, m2i = Mc, 
       sd1i = Se, sd2i = Sc, measure = "ROM")
m1 <- rma(yi = yi, sei = vi, method = "ML", 
          test = "knha", control=list(stepadj=0.5), data = d2)
