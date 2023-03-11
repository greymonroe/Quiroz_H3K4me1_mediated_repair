# Assuming no mutation bias (gene body and ingergenic region have the same mutation rate)
# x = the total number of mutations across genome
# d = the percentage decrease of ns mutations in gene body due to selection
# l = total length of genome (374471240)
# p = the percent of gene bodies that is coding sequences (0.3763199)
# r = expected ns/s ratio
# mu_g = the observed genic mutation rate
# mu_i = the observed intergenic mutation rate
# this is the observed muations in gene body:
#p*x*r/(r+1)*d + p*x*1/(r+1) = mu_g*p*l
# this is the observed mutations in intergenic region:
#(1-p)*x = mu_i*(1-p)*l
# based on the two equations:
#(r*d+1)/(r+1) = mu_g/mu_i
d = ((mu_g/mu_i)*(r+1)-1)/r
# if:
l=374471240
p=0.3763199
r <- 2.329456
mu_g<-9867/110152630 #genic mutaiton rate mutations per base pair
mu_i<-33616/(374471240-110152630) # non-genic

(33616/(374471240-110152630)-9867/110152630)/(33616/(374471240-110152630) )
# d = 0.664085
# meaning (1-d)% non-synonymous mutations in the gene body has been removed by selection
# assuming the expected ns/s is correct, what would be the observed ns and s ratio necessary to explain the observed reduction in gene body mutation rates:
(r/(r+1)*d)/(1/(r+1)) = 1.546957
# if the expected n/s is wrong, what would it have to be in order for the observed values to happen due to selection that reduces the gene body mutation rate as such
(r/(r+1)*d)/(1/(r+1)) = 5370/2155
(r*d+1)/(r+1) = mu_g/mu_i
# based on the two equations:
r = (5370/2155+1)/(mu_g/mu_i)-1
r = 7.348169

chisq.test(c(2155, 2155*1.546957), p=prop.table(c(1, 2.329456)))

