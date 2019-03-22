require(graphics)

pvalues_list <- list()
pvalues_list[[1]] <- c(0.001, 0.005, 0.1, 0.5, 0.01)
pvalues_list[[2]] <- c(0.03, 0.8, 0.47, 0.1)
pvalues_list[[3]] <- c(0.0003, 0.4, 0.002)

result = list()
for(i in 1:length(pvalues_list)){
    result[[i]] <- p.adjust(pvalues_list[[i]], "BH")
}

print(result)