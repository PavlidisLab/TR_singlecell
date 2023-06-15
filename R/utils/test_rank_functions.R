library(assertthat)


test_mat <- cor_l$Microglia[1:1000, 1:1000]


# Default is random
colrank1 <- colrank_mat(test_mat) 
colrank2 <- colrank_mat(test_mat) 
are_equal(colrank1, colrank2)


# Self as best rank=1
all(diag(colrank1) == 1)



# Default is random
allrank1 <- allrank_mat(test_mat)
allrank2 <- allrank_mat(test_mat)
are_equal(allrank1, allrank2)


