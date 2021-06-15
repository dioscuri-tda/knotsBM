setwd("~/GitHub/knotsBM_example/")

library(data.table) # to read csv faster

source('R/BallMapper.R')
sourceCpp('R/BallMapper.cpp')


#####################################################
#For standard, not symmetric data:
#This is the fucntion to create the standard BM
BallMapperCpp <- function( points , values , epsilon )
{
  output <- SimplifiedBallMapperCppInterface( points , values , epsilon )
  colnames(output$vertices) = c('id','size')
  return_list <- output
}#BallMapperCpp



# load the dataset
jones_13 <- fread('data/Jones_upto_13.csv', sep=',', header = TRUE)
sapply(jones_13[, 1:10], class)

# the table has 57 columns, the first 6 are info, from 7 to 57 there are the coefficients
colors <- jones_13[, 1:6]
coeff <- jones_13[, 7:57]
# add the norm of the coefficients
colors$norm <- wordspace::rowNorms(as.matrix(coeff))

# write the tables to file
write.table(colors, 'output/jones_13/Jones_upto_13_colors.csv')
write.table(coeff, 'output/jones_13/Jones_upto_13_coeff.csv')


# create a bm of radius epsilon, and color by the signature 
epsilon <- 20

print("COMPUTING BM")
start <- Sys.time()
jones_13_BM <- BallMapperCpp(coeff, jones_13$signature, epsilon)
print("DONE")
print(Sys.time() - start)
print("SAVING")
# save BM to file
# WARNING to not run more than once, it will corrupt the output files
storeBallMapperGraphInFile(jones_13_BM, filename = paste0("output/jones_13/", epsilon))
print("THE END")

# plot the BM
ColorIgraphPlot(jones_13_BM, seed=42)



#####################################################
# let's now consider the same dataset but this time we consider both a knot
# and its mirror
jones_13_MIRRORS <- fread('data/Jones_upto_13_MIRRORS.csv', sep=',', header = TRUE)
sapply(jones_13_MIRRORS[, 1:10], class)

# the table has 57 columns, the same structure as before, but twice the rows
colors <- jones_13_MIRRORS[, 1:6]
coeff <- jones_13_MIRRORS[, 7:57]
# add the norm of the coefficients
colors$norm <- wordspace::rowNorms(as.matrix(coeff))

# write the tables to file, warning do not run twice on the same output
write.table(colors, 'output/jones_13_MIRRORS/Jones_upto_13_MIRRORS_colors.csv')
write.table(coeff, 'output/jones_13_MIRRORS/Jones_upto_13_MIRRORS_coeff.csv')



# create a bm of radius epsilon, and color by the signature 
epsilon <- 20

print("COMPUTING BM")
start <- Sys.time()
jones_13_MIRRORS_BM <- BallMapperCpp(coeff, jones_13_MIRRORS$signature, epsilon)
print("DONE")
print(Sys.time() - start)
print("SAVING")
# save BM to file, warning do not run twice on the same output
storeBallMapperGraphInFile(jones_13_MIRRORS_BM, filename = paste0("output/jones_13_MIRRORS/", epsilon))
print("THE END")

# plot the BM, now there are also clusters with negative signature
ColorIgraphPlot(jones_13_MIRRORS_BM, seed=42)





#####################################################
# We know want to exploit the symmetry of the Jones polynomials with respect to the mirroring
# in the creation of the BM
# We use the following procedure

#For symmetric data:
BallMapperGroupActionCpp <- function( points , values , epsilon , orbit )
{
  output <- SimplifiedBallMapperCppInterfaceGroupAction( points , values , epsilon , orbit )
  colnames(output$vertices) = c('id','size')
  return_list <- output
}#BallMapperCpp

# we will again use jones_13_MIRRORS but we need also this file
# this files links each point (rows in jones_13_MIRRORS) with its mirror 
orbit <- read.csv('data/Jones_upto_13_MIRRORS_orbs.csv',header=FALSE)
View(orbit)

# we can then compute a Symmetric BM
# create a bm of radius epsilon, and color by the signature 
epsilon <- 20

print("COMPUTING SYMMETRIC BM")
start <- Sys.time()
jones_13_MIRRORS_SYMBM <- BallMapperGroupActionCpp(coeff, jones_13_MIRRORS$signature, epsilon, orbit)
print("DONE")
print(Sys.time() - start)
print("SAVING")
# save BM to file
storeBallMapperGraphInFile(jones_13_MIRRORS_SYMBM, filename = paste0("output/jones_13_MIRRORS/SYM_", epsilon))
print("THE END")

# plot the BM, now flares with opposite signature are exactly the same
ColorIgraphPlot(jones_13_MIRRORS_SYMBM, seed=42)







#####################################################
# We know want to exploit the symmetry of the Jones polynomials with respect to the mirroring
# in the creation of the BM
# We use the following procedure
# JONES upto13n from Khovanov data (for mapper on BM)

#For symmetric data:
BallMapperGroupActionCpp <- function( points , values , epsilon , orbit )
{
  output <- SimplifiedBallMapperCppInterfaceGroupAction( points , values , epsilon , orbit )
  colnames(output$vertices) = c('id','size')
  return_list <- output
}#BallMapperCpp


jones_13n_MIRRORS <- fread('data/Jones_fromK_upto_13n_MIRRORS.csv', sep=',', header = TRUE)
sapply(jones_13n_MIRRORS[, 1:10], class)

# the table has 57 columns, the same structure as before, but twice the rows
colors <- jones_13n_MIRRORS[, 1:5]
coeff <- jones_13n_MIRRORS[, 6:40]
# add the norm of the coefficients
colors$norm <- wordspace::rowNorms(as.matrix(coeff))

# write the tables to file, warning do not run twice on the same output
write.table(colors, 'output/jones_fromK_upto_13n_MIRRORS/Jones_fromK_upto_13n_MIRRORS_colors.csv')
write.table(coeff, 'output/jones_fromK_upto_13n_MIRRORS/Jones_fromK_upto_13n_MIRRORS_coeff.csv')



# we will again use jones_13_MIRRORS but we need also this file
# this files links each point (rows in jones_13_MIRRORS) with its mirror 
orbit <- read.csv('data/Jones_fromK_upto_13n_MIRRORS_orbs.csv',header=FALSE)
View(orbit)

# we can then compute a Symmetric BM
# create a bm of radius epsilon, and color by the signature 
epsilon <- 20

print("COMPUTING SYMMETRIC BM")
start <- Sys.time()
jones_13n_MIRRORS_SYMBM <- BallMapperGroupActionCpp(coeff, jones_13n_MIRRORS$signature, epsilon, orbit)
print("DONE")
print(Sys.time() - start)
print("SAVING")
# save BM to file
storeBallMapperGraphInFile(jones_13n_MIRRORS_SYMBM, filename = paste0("output/jones_fromK_upto_13n_MIRRORS/SYM_", epsilon))
print("THE END")

# plot the BM, now flares with opposite signature are exactly the same
ColorIgraphPlot(jones_13n_MIRRORS_SYMBM, seed=42)


khov_13n_MIRRORS <- fread('data/Khovanov_upto_13n_MIRRORS.csv', sep=',', header = TRUE)
sapply(khov_13n_MIRRORS[, 1:10], class)
colors_K <- khov_13n_MIRRORS[, 1:5]
coeff_K <- khov_13n_MIRRORS[, 6:1922]
write.table(coeff_K, 'output/jones_fromK_upto_13n_MIRRORS/Khovanov_upto_13n_MIRRORS_coeff.csv')
write.table(colors_K, 'output/jones_fromK_upto_13n_MIRRORS/Khovanov_upto_13n_MIRRORS_colors.csv')


