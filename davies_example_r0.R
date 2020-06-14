# Davies example R0
setwd("~/Vaccine Strategy")

matrices = readRDS("all_matrices.rds")
names(matrices) # to see your options

mat = matrices[["United States of America"]]

cm = mat$home + mat$work + mat$school + mat$other

u     <- c(0.33, 0.33, 0.37, 0.37, 0.69, 0.69, 0.81, 0.81, 0.74, 0.74, 0.8, 0.8, 0.89, 0.89, 0.77, 0.77) # Ref: Davies
u <- u/100

Du <- diag(u)
Dy <- diag(5, 16)

#ngm <- u * cm * 5

ngm <- Du %*% cm %*% Dy
R0 <- abs(eigen(ngm)$values[1])
