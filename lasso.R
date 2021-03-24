source('clean_data.R')
source('boot_glmnet.R')
library(Matrix)

d0 <- d %>% reshape2::melt(id.vars = c("Batch",
                                       "Time",
                                       "Sex",
                                       "PG.Genes",
                                       "tbi")) %>%
  dplyr::mutate(animal_id = factor(PG.Genes),
                Sex = factor(Sex)) %>%
  dplyr::rename(protein = variable)

f <- as.formula(value ~ Time*protein)
y <- d0$value
# Second step: using model.matrix to take advantage of f
x <- sparse.model.matrix(f, d0)[, -1]

R <- 100  # n bootstrap replications - why so many? Lots of data
best_alpha <- 1 
best_lambda <- .004

m1 <- boot_glmnet(x = x,
                  y = y,
                  alpha = 1,
                  lambda = .004, 
                  R = 25000, verbose = TRUE)
saveRDS(m1, 'boot_m1.rds')
m2 <- boot_glmnet(x = x,
                 y = y,
                 alpha = 1,
                 lambda = .004, 
                 R = 25000, verbose = TRUE)
saveRDS(m2, 'boot_m2.rds')
boot_combine(m1, m2)

m3 <- quant.boot_glmnet(m3, ci = .95)
plot(m3, p.lasso.thresh = .95)
summary(m3, p.lasso.thresh = .95)



# p.lasso.thresh is the minimum threshold 
# for the proportion of times the LASSO operator
# estimated a non-zero coefficient for a variable
# during bootstrap replications
# these might have CIs which overlap zero but they
# are nevertheless quite likely to be important variables

summary(m1, p.lasso.thresh = .95)
plot_coefs(m1, p.lasso.thresh = 1)

saveRDS(m1, file = 'm1_boot_parallel.rds')
