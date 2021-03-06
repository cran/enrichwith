#' Create a enrichwith skeleton
#'
#'
#' Create an enrichwith skeleton file for the structured
#' implementation of methods to compute new components for objects of
#' a specific class
#'
#' @param class the class of the objects to be enriched
#' @param option a character vector with components the enrichment
#'     options
#' @param description a character vector of length
#'     \code{length(options)} with components the description of the
#'     enrichment options
#' @param component a list of as many character vectors as
#'     \code{length(option)}, specifying the names of the components
#'     that each option will add to the object after enrichment
#' @param path the path where the skeleton file will be created
#' @param filename the name of the skeleton file
#' @param attempt_rename attempt to rename syntactically incorrect
#'     component names? Default is \code{TRUE}
#'
#' @return A file with the necessary functions to use
#'     \code{enrichwith} infrastructure. The skeleton consists of the
#'     following functions
#' \itemize{
#'
#' \item One \code{compute_component.class} function per component
#' name from \code{unique(unlist(component))}. The function takes as
#' input the object to be enriched and returns as output the component
#' to be added to the object.
#'
#' \item The \code{get_enrichment_options.class} function, that takes
#' as input the object to be enriched and an enrichment option, and
#' returns the names of the components that will be appended to the
#' object for this option. This function can also be used to list the
#' available options and print their description.
#'
#' \item The \code{enrich.class} function
#'
#' }
#' @examples
#' \dontrun{
#' # Set the directory where the skeleton is placed
#' my_path <- "~/Downloads"
#' # This is the call that created the enrichment skeleton for glms
#' # that ships with the package
#' create_enrichwith_skeleton(class = "glm",
#'       option = c("auxiliary functions", "score vector",
#'        "mle of dispersion", "expected information",
#'        "observed information", "first-order bias"),
#'       description = c("various likelihood-based quantities
#'         (gradient of the log-likelihood, expected and observed
#'         information matrix and first term in the expansion of
#'         the bias of the mle) and a simulate method as functions
#'         of the model parameters",
#'        "gradient of the log-likelihood at the mle",
#'        "mle of the dispersion parameter",
#'        "expected information matrix evaluated at the mle",
#'        "observed information matrix evaluated at the mle",
#'        "first term in the expansion of the bias of the mle
#'         at the mle"),
#'       component = list("auxiliary_functions", "score_mle",
#'                        "dispersion_mle",
#'                        "expected_information_mle",
#'                        "observed_information_mle",
#'                        "bias_mle"),
#'       path = my_path,
#'       attempt_rename = FALSE)
#'
#' }
#' @export
create_enrichwith_skeleton <- function(class,
                                      option,
                                      description,
                                      component,
                                      path,
                                      filename = paste0(class, "_options.R"),
                                      attempt_rename = TRUE) {

    con <- file(paste(path, filename, sep = "/"), open = "w")

    o_length <- length(option)
    d_length <- length(description)

    message("* Checking supplied class for validity.")
    class_input <- as.character(class)
    ## Just making sure that a syntactically correct name is used for the class
    class_valid <- make.names(class, unique = TRUE, allow_ = FALSE)

    message("* Checking supplied class for validity.")
    option <- as.character(option)

    message("* Checking supplied descriptions for validity.")
    description <- as.character(description)
    if (d_length != o_length) {
        stop("ength(description) is not equal to length(option)")
    }

    message("* Checking supplied component for validity.")
    if (length(component) != o_length) {
        stop("length(component) is not equal to length(option)")
    }
    ## So even if there is only one component, component will become a list
    component <- lapply(component, function(comp) {
        comp_input <- as.character(comp)
        comp_valid <- make.names(comp, unique = TRUE, allow_ = FALSE)
        valid <- comp_input == comp_valid
        if (!all(valid)) {
            if (attempt_rename) {
                warning(gettextf("not syntactically valid component names: %s were renamed to %s",
                                 paste0(paste0("'", comp_input[!valid], "'"), collapse = ", "),
                                 paste0(paste0("'", comp_valid[!valid], collapse = ", "))))
            }
            else {
                ## close(con)
                warning(gettextf("not syntactically valid component names: %s",
                                 paste0(paste0("'", comp_input[!valid], "'"), collapse = ", ")))
            }
        }
        if (attempt_rename) comp_valid else comp_input
    })
    dat <- list(class = class_input,
                option = gsub("\"","'", deparse(option, width.cutoff = 500)),
                description = gsub("\"","'", deparse(description, width.cutoff = 500)),
                component = gsub("\"","'", deparse(component, width.cutoff = 500)))

    ## enrich function
    message("* Setting up enrich.", class_input, " function")
    ## Get _options template
    template_path <- system.file("templates", "template_enrich.R", package = "enrichwith",  mustWork = TRUE)
    ## Read template and replace according to class/options/description/component
    template_out <- whisker::whisker.render(readLines(template_path), data = dat)
    ## Write into _options file
    writeLines(template_out, con = con)

    ## get_enrichment_options fuction
    message("* Setting up get_enrichment_options.", class_input, " function")
    ## Get _options template
    template_path <- system.file("templates", "template_get_enrichment_options.R", package = "enrichwith",  mustWork = TRUE)
    ## Read template and replace according to class/options/description/component
    template_out <- whisker::whisker.render(readLines(template_path), data = dat)
    ## Write into _options file
    writeLines(template_out, con = con)

    ## compute_component functions
    message("* Setting up compute_component.", class_input, " function")
    components <- unique(unlist(component))
    template_path <- system.file("templates", "template_compute_component.R", package = "enrichwith",  mustWork = TRUE)
    for (j in seq.int(length(components))) {
        template_out <-  whisker::whisker.render(readLines(template_path),
                                                 data = list(component = components[j],
                                                             class = class_input))
        writeLines(template_out, con = con)
    }

    call <- match.call()

    writeLines(paste("\n\n\n\n", "## ## Call that produced the enrichwith template for the current script:\n", paste("##", deparse(call), collapse = "\n"), "\n", collapse = "", sep = ""), con = con)


    close(con)


}


