# Get IFR

# meta-regression from Levin et al.
# https://www.medrxiv.org/content/10.1101/2020.07.23.20160895v4.full.pdf+html

age = seq(0,89,1)
IFR_byones = exp(-7.56 + 0.121*age)
IFR_bytens = rep(0, 9)

count = 1  
for (i in 1:9){
  IFR_bytens[i] = sum(IFR_byones[count:(count+9)])/10
  count = count + 10
}