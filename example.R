# Download COBRE data from "graphclass"
#library(devtools)
#devtools::install_github("jesusdaniel/graphclass")
library(graphclass)
library(docstring)

data("COBRE.data")

X = COBRE.data$X.cobre
y = COBRE.data$Y.cobre

# convert the data to adjacency matrices
m = nrow(X)
p = ncol(X)
n = 263 # number of nodes in the COBRE data

A_list = lapply(1:m, function(i) get_matrix(X[i, ]))


# Power parcelation communities
data("power.parcellation")
power.parcellation$Master.Assignment
# Create membership matrix
Z_power <- Z_membership_from_vector(power.parcellation$Master.Assignment)
Z_power = Z_power[-75,] # delete missing node
# print community sizes
colSums(Z_power)

# plot a network with Power atlas
plot_adjmatrix(A_list[[1]], 
               communities = Z_membership_to_list(Z_power),
               community_labels = c(1:13, -1), axislabel = "brain systems")

# load supervised community detection function -----------------
source("glmblock/load_glmblock.R")
docstring(glmblock)
docstring(glmblock.fit)


# 1) Fit solution of Power parcellation --------------------------------------------
# --fit constrained regression with Power communities
glm.power <- glmblock.fit(Adj_list = A_list, y = y,
                          family = "binomial", intercept = TRUE,
                          Z = Z_power, lambda = 0, gamma = 1e-2)
# --plot solution
plot_adjmatrix(glm.power$B, communities = Z_membership_to_list(Z_power),
               community_labels = c(1:13, -1), axislabel = "brain systems")




# 2) Find spectral clustering solution -------------------------------------------

glm.SCD.SC = glmblock(Adj_list = A_list, y = y, 
                      family = "binomial", intercept = TRUE,
                      K = 14, gamma = 1e-2, lambda = 0, 
                      method.communities = "SC")


# --plot solution with Power parcelation
plot_adjmatrix(glm.SCD.SC$glm.result$B, 
               communities = Z_membership_to_list(Z_power),
               community_labels = c(1:13, -1), axislabel = "brain systems")




# 3) Find ADMM solution -----------------------------------------------------------

glm.SCD.ADMM = glmblock(Adj_list = A_list, y = y, 
                      family = "binomial", intercept = TRUE,      # regression parameters
                      K = 14, gamma = 1e-2, lambda = 0,        # penalty parameters
                      method.communities = "ADMM",  # ADMM parameters
                      verbose = TRUE)

plot_adjmatrix(glm.SCD.ADMM$glm.result$B, 
               communities = Z_membership_to_list(Z_power),
               community_labels = c(1:13, -1), axislabel = "brain systems")


# Compare fit of different solutions, smallest values of the objective function are better fit to the training data
glm.power$objective
glm.SCD.SC$glm.result$objective
glm.SCD.ADMM$glm.result$objective
