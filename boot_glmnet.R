bootfit <- function(x, y, alpha, lambda, t0, R, z, verbose, ...){
  require(glmnet)
  if(verbose){cat('\rFitting bootstrap replicate',
                  z,
                  'of',
                  R)
  }
  k <- sample(1:length(y), length(y), replace = TRUE)
  y <- y[k]
  x <- x[k,]
  m <- glmnet(x = x, y = y, alpha = alpha, lambda = lambda, family = "gaussian")
  pred <- predict(m, x, lambda = lambda, type = "response")
  t1 <- Sys.time()
  elapsed <- as.numeric(t1-t0)
  tms <- ifelse(elapsed/z*R < 120, 1, 60)
  tms.label <- ifelse(elapsed/z*R < 120, 'sec', 'min')
  if(verbose){cat('\b . Time elapsed',
                  round(elapsed/tms, 1), tms.label,
                  '\b. Expected time to complete', 
                  round((elapsed/z)/tms*R, 1), tms.label, '\b.             ')
  }
  return(list(RMSE = caret::RMSE(pred, y),
              R2 = caret::R2(pred, y)[1,1],
              coefs = coef(m, lambda = lambda)))
}


boot_glmnet <- function(x, y, alpha, lambda, verbose = TRUE, R = 100, ...){
  t0 <- Sys.time()
  if(verbose) {cat('\n')}
  out <- lapply(1:R, function(z) bootfit(y = y,
                                                x = x,
                                                alpha = alpha,
                                                lambda = lambda,
                                                t0 = t0,
                                                R = R,
                                                z = z,
                                                verbose = verbose,
                                                ...))
  out <- simplify2array(out)
  
  r2vec <- unlist(out['R2',])
  rmsevec <- unlist(out['RMSE',])
  coefs <- out['coefs',]
  coefmat <- do.call(cbind, coefs)
  obj <- list(alpha = alpha,
              lambda = lambda,
              R = R,
              r2vec = r2vec,
              rmsevec = rmsevec,
              coefmat = coefmat)
  class(obj) <- 'boot_glmnet'
  return(obj)
}

boot_combine <- function(fit1, fit2){
  if(fit1$alpha != fit2$alpha){
    warning('\nModel alphas not equal, cannot combine')
    return(NULL)
  }
  if(fit1$lambda != fit2$lambda){
    warning('\nModel lambdas not equal, cannot combine')
    return(NULL)
  }
  if(nrow(fit1$coefmat) != nrow(fit2$coefmat)){
    warning('\nModel matrixes not equal, cannot combine')
    return(NULL)
  }
  r2vec <- c(fit1$r2vec, fit2$r2vec)
  rmsevec <- c(fit1$rmsevec, fit2$rmsevec)
  coefmat <- cbind(fit1$coefmat, fit2$coefmat)
  obj <- list(alpha = fit1$alpha,
              lambda = fit1$lambda,
              R = fit1$R + fit2$R,
              r2vec = r2vec,
              rmsevec = rmsevec,
              coefmat = coefmat)
  class(obj) <- 'boot_glmnet'
  return(obj)
  
}

quant.boot_glmnet <- function(obj, ci = .95){
  if(ci > 1){
    ci <- ci/100
    warning('\nci > 1, assumed to be percentage (divided by 100)\n')
  }
  ci.lwr <- (1-abs(ci))/2
  ci.upr <- 1-((1-abs(ci))/2)
  R <- obj$R
  coefmat <- obj$coefmat
  r2vec <- obj$r2vec
  rmsevec <- obj$rmsevec
  
  coefs.lwr <- apply(coefmat, 1, function(z) quantile(z, ci.lwr, na.rm = TRUE))
  coefs.upr <- apply(coefmat, 1, function(z) quantile(z, ci.upr, na.rm = TRUE))
  coefs.m <- apply(coefmat, 1, function(z) mean(z, na.rm = TRUE))
  # p(selection|alpha, lambda)
  p.lasso <- apply(coefmat, 1, function(z) sum(abs(z)>0)/R)
  # Approx CIs for p.lasso given R replicates
  p.lasso.lwr <- p.lasso - qnorm(ci.upr)*sqrt(p.lasso*(1-p.lasso)/R)
  p.lasso.upr <- p.lasso + qnorm(ci.upr)*sqrt(p.lasso*(1-p.lasso)/R)
  p.lasso.upr[p.lasso.upr > 1] <- 1
  p.lasso.lwr[p.lasso.lwr < 0] <- 0
  r2 <- c(mean(r2vec, na.rm = TRUE), quantile(r2vec, c(ci.lwr, ci.upr), na.rm = TRUE))
  rmse <- c(mean(rmsevec, na.rm = TRUE), quantile(rmsevec, c(ci.lwr, ci.upr), na.rm = TRUE))
  coefs <- data.frame(coefs.m, coefs.lwr, coefs.upr,
                      p.lasso, p.lasso.lwr, p.lasso.upr)
  obj <- list(r2 = r2,
              rmse = rmse,
              coefs = coefs,
              alpha = obj$alpha,
              lambda = obj$lambda,
              R = R,
              ci = ci)
  class(obj) <- 'boot_glmnet_q'
  return(obj)
}

# extract data.frame of P(selection|alpha, lambda) and a CI
p.lasso <- function(obj, ...){
  if(!class(obj) %in% c('boot_glmnet', 'boot_glmnet_q')){
    warning('\nNot a boot_glmnet or boot_glmnet_q object')
  }
  x <- obj
  if(class(obj) == 'boot_glmnet'){
    x <- quant.boot_glmnet(obj, ...)
  }
  
  y <- x$coefs[, c(4:6)]
  colnames(y) <- c('p.lasso', paste0('lwr.', x$ci*100), paste0('upr.', x$ci*100))
  attr(y, 'ci') <- x$ci
  attr(y, "R") <- x$R
  attr(y, 'alpha') <- x$alpha
  attr(y, 'lambda') <- x$lambda
  return(y)
}

coef.boot_glmnet <- function(obj, ...){
  x <- quant.boot_glmnet(obj, ...)
  y <- x$coefs[, c(1:3)]
  colnames(y) <- c('Estimate', paste0('lwr.', x$ci*100), paste0('upr.', x$ci*100))
  attr(y, 'ci') <- x$ci
  attr(y, "R") <- x$R
  attr(y, 'alpha') <- x$alpha
  attr(y, 'lambda') <- x$lambda
  return(y)
}

coef.boot_glmnet_q <- function(obj, ...){
  x <- obj
  y <- x$coefs[, c(1:3)]
  colnames(y) <- c('Estimate', paste0('lwr.', x$ci*100), paste0('upr.', x$ci*100))
  attr(y, 'ci') <- x$ci
  attr(y, "R") <- x$R
  attr(y, 'alpha') <- x$alpha
  attr(y, 'lambda') <- x$lambda
  return(y)
}

sig_coefs <- function(obj, p.lasso.thresh = 1, ...){
  if(class(obj) == 'boot_glmnet'){
    obj <- quant.boot_glmnet(obj, ...)
  }
  x <- obj$coefs
  flag <- sign(x[,2]) == sign(x[,3])
  flag[x[,2] == 0] <- FALSE
  plt <- x$p.lasso > p.lasso.thresh
  
  names(x) <- c('Estimate', 
                 paste0('lwr'),
                 paste0('upr'),
                 'p.lasso',
                 paste0('p.lasso.lwr', obj$ci*100),
                 paste0('p.lasso.upr', obj$ci*100))
  return(data.frame(x, flag = flag, p.lasso.thresh = plt))
}


summary.boot_glmnet <- function(obj, p.lasso.thresh = 1, digits = 3, ...){
  obj <- quant.boot_glmnet(obj, ...)
  cat('\nalpha: ')
  print(obj$alpha)
  cat('\nlambda: ')
  print(obj$lambda)
  
  
  cat('\nRMSE: \n')
  print(obj$rmse, digits = digits)
  cat('\nR-squared:  \n')
  print(obj$r2, digits = digits)
  ci <- sig_coefs(obj, p.lasso.thresh)
  cat('\nCoefficients with significant coefficients at p <', 1-obj$ci,'or P(selection with LASSO) >', p.lasso.thresh, ':  \n')
  ci$sig <- ifelse(ci$flag, '*', ' ')
  print(ci[ci$flag | ci$p.lasso.thresh,-c(7,8)], digits = digits)
  cat('\n\n')
}

summary.boot_glmnet_q <- function(obj, p.lasso.thresh = 1, digits = 3, ...){
  cat('\nalpha: ')
  print(obj$alpha)
  cat('\nlambda: ')
  print(obj$lambda)
  
  
  cat('\nRMSE: \n')
  print(obj$rmse, digits = digits)
  cat('\nR-squared:  \n')
  print(obj$r2, digits = digits)
  ci <- sig_coefs(obj, p.lasso.thresh)
  cat('\nCoefficients with significant coefficients at p <', 1-obj$ci,'or P(selection with LASSO) >', p.lasso.thresh, ':  \n')
  ci$sig <- ifelse(ci$flag, '*', ' ')
  print(ci[ci$flag | ci$p.lasso.thresh,-c(7,8)], digits = digits)
  cat('\n\n')
}

plot.boot_glmnet <- function(obj, keep = NULL, ...){
  s <- sig_coefs(obj, ...) %>% tibble::rownames_to_column() %>%
    rename(coef = rowname) %>%
    filter(flag | coef %in% keep | p.lasso.thresh, coef != '(Intercept)') %>%
    mutate(coef = forcats::fct_reorder(coef, Estimate))
  
  g <- ggplot(s, aes(x = coef, y = Estimate, colour = flag)) +
    geom_errorbar(width = .2, aes(ymin = lwr, ymax = upr)) +
    geom_point() +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    scale_colour_manual(values = c("black", "red")) +
    xlab("") +
    ylab("[Penalized] coefficient (in standardized units, SD)") +
    coord_flip() +
    theme_bw() +
    theme(legend.position = "none")
  g
}

plot.boot_glmnet_q <- function(obj, keep = NULL, ...){
  s <- sig_coefs(obj) %>% tibble::rownames_to_column() %>%
    rename(coef = rowname) %>%
    filter(flag | coef %in% keep | p.lasso.thresh, coef != '(Intercept)') %>%
    mutate(coef = forcats::fct_reorder(coef, Estimate))
  
  g <- ggplot(s, aes(x = coef, y = Estimate, colour = flag)) +
    geom_errorbar(width = .2, aes(ymin = lwr, ymax = upr)) +
    geom_point() +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    scale_colour_manual(values = c("black", "red")) +
    xlab("") +
    ylab("[Penalized] coefficient (in standardized units, SD)") +
    coord_flip() +
    theme_bw() +
    theme(legend.position = "none")
  g
}

mc_boot_glmnet <- function(blocks = 4L,
                           x,
                           y,
                           alpha = 1,
                           lambda,
                           R = 100,
                           verbose = TRUE){
  if(blocks < 2){
    warning('\nblocks must be >=2, increasing blocks to 2\n')
    blocks = 2L
  }
  require(future.apply)
  plan(multisession)
  k <- future_replicate(blocks,
              boot_glmnet(x = x,
                          y = y,
                          alpha = alpha,
                          lambda = lambda,
                          R = R, 
                          verbose = verbose,
                          ...),
              simplify = "list",
              future.seed = TRUE)
  j <- k[,1]
  for(i in 2:blocks){
    j <- boot_combine(j, k[,i])
  }
  return(j)
}
