#####################################################
# MAPPER ON BALL MAPPER
# TO USE BEFORE mapper_on_BM.ipynb
#####################################################
# JONES upto13n from Khovanov data (for mapper on BM)

setwd("~/GitHub/knotsBM_example/")

library(data.table) # to read csv faster

library(Rcpp)
sourceCpp('R/BallMapper.cpp')

source('R/BallMapper_utils.R')


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