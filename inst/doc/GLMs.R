## ---- echo = TRUE, eval = TRUE-------------------------------------------
clotting <- data.frame(conc = c(5,10,15,20,30,40,60,80,100,
                                5,10,15,20,30,40,60,80,100),
                       time = c(118, 58, 42, 35, 27, 25, 21, 19, 18,
                                69, 35, 26, 21, 18, 16, 13, 12, 12),
                       lot = factor(c(rep(1, 9), rep(2, 9))))

## ---- echo = TRUE, eval = TRUE-------------------------------------------
library("ggplot2")
clottingML <- glm(time ~ log(conc) * lot, family = Gamma, data = clotting)
alpha <- 0.01
pr_out <- predict(clottingML, type = "response", se.fit = TRUE)
new_data <- clotting
new_data$time <- pr_out$fit
new_data$type <- "fitted"
clotting$type <- "observed"
all_data <- rbind(clotting, new_data)
new_data <- within(new_data, {
    low <- pr_out$fit - qnorm(1 - alpha/2) * pr_out$se.fit
    upp <- pr_out$fit + qnorm(1 - alpha/2) * pr_out$se.fit
})
ggplot(all_data) + geom_point(aes(conc, time, col = type), alpha = 0.8) +
    geom_segment(data = new_data, aes(x = conc, y = low, xend = conc, yend = upp, col = type)) +
    facet_grid(. ~ lot) + theme_bw() + theme(legend.position = "top")

## ---- echo = TRUE, eval = TRUE-------------------------------------------
library("enrichwith")
enriched_clottingML <- enrich(clottingML, with = "auxiliary functions")
scores_clottingML <- enriched_clottingML$auxiliary_functions$score

## ---- echo = TRUE, eval = TRUE-------------------------------------------
scores_clottingML <- get_score_function(clottingML)

## ---- echo = TRUE, eval = TRUE-------------------------------------------
scores_clottingML()

## ---- echo = TRUE, eval = TRUE-------------------------------------------
info_clottingML <- enriched_clottingML$auxiliary_functions$information

## ---- echo = TRUE, eval = TRUE-------------------------------------------
info_clottingML <- get_information_function(clottingML)

## ---- echo = TRUE, eval = TRUE-------------------------------------------
summary_clottingML <- summary(clottingML)

## ---- echo = TRUE, eval = TRUE-------------------------------------------
summary_std_errors <- coef(summary_clottingML)[, "Std. Error"]
einfo <- info_clottingML(dispersion = summary_clottingML$dispersion)
all.equal(sqrt(diag(solve(einfo)))[1:4], summary_std_errors, tolerance = 1e-05)

## ---- echo = TRUE, eval = TRUE-------------------------------------------
oinfo <- info_clottingML(dispersion = summary_clottingML$dispersion, type = "observed")
all.equal(oinfo[1:4, 1:4], einfo[1:4, 1:4])

## ---- echo = TRUE, eval = TRUE-------------------------------------------
clottingML_log <- update(clottingML, family = Gamma("log"))
summary_clottingML_log <- summary(clottingML_log)
info_clottingML_log <- get_information_function(clottingML_log)
einfo_log <- info_clottingML_log(dispersion = summary_clottingML_log$dispersion, type = "expected")
oinfo_log <- info_clottingML_log(dispersion = summary_clottingML_log$dispersion, type = "observed")
round(einfo_log, 3)
round(oinfo_log, 3)

## ---- echo = TRUE, eval = TRUE-------------------------------------------
clottingML_nested <- update(clottingML, . ~ log(conc))
enriched_clottingML_nested <- enrich(clottingML_nested, with = "mle of dispersion")
coef_full <- coef(clottingML)
coef_hypothesis <- structure(rep(0, length(coef_full)), names = names(coef_full))
coef_hypothesis_short <- coef(enriched_clottingML_nested, model = "mean")
coef_hypothesis[names(coef_hypothesis_short)] <- coef_hypothesis_short
disp_hypothesis <- coef(enriched_clottingML_nested, model = "dispersion")
scores <- scores_clottingML(coef_hypothesis, disp_hypothesis)
info <- info_clottingML(coef_hypothesis, disp_hypothesis)

## ---- echo = TRUE, eval = TRUE-------------------------------------------
(score_statistic <- drop(scores%*%solve(info)%*%scores))

## ---- echo = TRUE, eval = TRUE-------------------------------------------
pchisq(score_statistic, 2, lower.tail = FALSE)

## ---- echo = TRUE, eval = TRUE-------------------------------------------
coef_full[3:4]%*%solve(solve(info)[3:4, 3:4])%*%coef_full[3:4]

## ---- echo = TRUE, eval = TRUE-------------------------------------------
(deviance(clottingML_nested) - deviance(clottingML))/disp_hypothesis

## ---- echo = TRUE, eval = TRUE-------------------------------------------
simulate_clottingML <- get_simulate_function(clottingML)
simulate_clottingML(nsim = 3, seed = 123)

## ---- echo = TRUE, eval = TRUE-------------------------------------------
simulate(clottingML, nsim = 3, seed = 123)

## ---- echo = TRUE, eval = TRUE-------------------------------------------
coefficients <- c(0, 0.01, 0, 0.01)
dispersion <- 0.001
samples <- simulate_clottingML(coefficients = coefficients, dispersion = dispersion, nsim = 100000, seed = 123)

## ---- echo = TRUE, eval = TRUE-------------------------------------------
means <- 1/(model.matrix(clottingML) %*% coefficients)
variances <- dispersion * means^2
max(abs(rowMeans(samples) - means))
max(abs(apply(samples, 1, var) - variances))

## ---- echo = TRUE, eval = TRUE, cache = TRUE-----------------------------
compute_pvalues <- function(n, parameter, coefficients, dispersion, nsimu = 100, verbose = FALSE) {
    require("plyr")
    require("doMC")
    ## Concentration grid
    conc <- seq(5, 100, length.out = n)
    ## A data frame that sets the design. What the response is does
    ## not matter here
    clotting_temp <- data.frame(conc = rep(conc, 2),
                                lot = factor(c(rep(1, n), rep(2, n))),
                                time = rgamma(2 * n, 2, 2))
    ## Fit a dummy GLM and then get a simulate method out of it
    clotting_temp_ML<- glm(time ~ log(conc) * lot, family = Gamma, data = clotting_temp)
    simulate_clotting <- get_simulate_function(clotting_temp_ML)
    pvalues_out <- NULL
    for (which in seq.int(length(parameter))) {
        if (verbose) {
            cat("setting", which, "out of", length(parameter), "\n")
        }
        parameter_setting <- coefficients
        parameter_setting[4] <- parameter[which]
        simu_samples <- simulate_clotting(parameter_setting,
                                          dispersion = dispersion,
                                          nsim = nsimu,
                                          seed = 123)
        results <- adply(simu_samples, 2, function(response) {
        clotting_temp$response <- unlist(response)
        ## Fit the full model for the current response vector and enrich it
        cfit <-  enrich(update(clotting_temp_ML, response ~ .), with = "auxiliary functions")
        ## Fit the nested model
        cfit_nested <- update(cfit, . ~ log(conc) + lot)
        enriched_cfit_nested <- enrich(cfit_nested, with = "mle of dispersion")
        coef_full <- coef(cfit)
        ## Prepare the vectors of the constrained estimates
        coef_hypothesis <- structure(rep(0, length(coef_full)), names = names(coef_full))
        coef_hypothesis_short <- coef(enriched_cfit_nested, model = "mean")
        coef_hypothesis[names(coef_hypothesis_short)] <- coef_hypothesis_short
        disp_hypothesis <- coef(enriched_cfit_nested, model = "dispersion")
        ## Compute score and information of the full mode at the constrained estimates
        scores <- get_score_function(cfit)(coef_hypothesis, disp_hypothesis)
        info <- get_information_function(cfit)(coef_hypothesis, disp_hypothesis)
        ## Compute statistics
        score_statistic <- drop(scores%*%solve(info)%*%scores)
        lr_statistic <- (deviance(cfit_nested) - deviance(cfit))/disp_hypothesis
        wald_statistic <- coef_full[4]%*%solve(solve(info)[4, 4])%*%coef_full[4]
        ## power
        data.frame(pvalues = c(pchisq(score_statistic, 1, lower.tail = FALSE),
                               pchisq(lr_statistic, 1, lower.tail = FALSE),
                               pchisq(wald_statistic, 1, lower.tail = FALSE)),
                   statistics = c(score_statistic, lr_statistic, wald_statistic),
                   method = factor(c("score", "lr", "wald")))
        }, .parallel = TRUE)
        pvalues_out <- rbind(pvalues_out, data.frame(results, parameter = parameter[which], sample_size = n))
    }
    pvalues_out
}

## ---- echo = TRUE, eval = FALSE, cache = TRUE----------------------------
#  library("dplyr")
#  registerDoMC(2)
#  
#  nsimu <- 500
#  coefficients <- coef(enriched_clottingML, "mean")
#  dispersion <- coef(enriched_clottingML, "dispersion")
#  parameters <- seq(0, 0.002, length = 20)
#  
#  ## Compute pvalues for each resolution of percentage concentration
#  pvalues <- NULL
#  for (n in c(9, 18, 36, 72)) {
#      set.seed(123)
#      pvalues <- rbind(pvalues, compute_pvalues(n, parameters, coefficients, dispersion, nsimu))
#  }
#  
#  ## Compute power for each combination of parameter, method and sample size
#  power_values_1 <- data.frame(pvalues %>%
#                               group_by(parameter, method, sample_size) %>%
#                               summarize(power = mean(pvalues < 0.01)),
#                               alpha = 0.01)
#  power_values_5 <- data.frame(pvalues %>%
#                               group_by(parameter, method, sample_size) %>%
#                               summarize(power = mean(pvalues < 0.05)),
#                               alpha = 0.05)
#  power_values <- rbind(power_values_1, power_values_5)
#  
#  ## Plot the power curves
#  ggplot(power_values) + geom_line(aes(x = parameter, y = power, col = sample_size, group = sample_size)) + facet_grid(method ~ alpha)

