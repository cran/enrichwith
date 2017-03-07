
glmplus <- function(...) {
    fit <- glm(...)
    fit <- enrich(fit, with = "all")

}
