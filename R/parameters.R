parameter <- function(name,
                      type = c("integer", "float", "logical", "mc", "file"),
                      is_null_possible = FALSE,
                      positive = FALSE,
                      nonnegative = FALSE,
                      mc = NULL,
                      file = NULL) {
    structure(list(name = name,
                   is_null_possible = is_null_possible,
                   positive = positive,
                   nonnegative = nonnegative,
                   mc = mc,
                   file = file),
              class = c(paste0(type, "_parameter"), parameter_class))
}

check_parameter_class <- function(parameter) {
    if (!inherits(parameter, parameter_class)) {
        stop("Invalid argument")
    }
}

check_parameter <- function(param, value) {
    if (is.null(value)) {
        if (param$is_null_possible) {
            return()
        } else {
            stop(paste("Parameter", param$name, "mustn't be null"))
        }
    }
    if (length(value) == 0) {
        stop(paste("Parameter", param$name, "mustn't be a zero length vector"))
    }
    UseMethod("check_parameter")
}

check_parameter.integer_parameter <- function (param, value) {
    value <- as.integer(value)
    if (is.na(value)) {
        stop(paste("Parameter", param$name, "mustn't be NA"))
    }
    if (param$positive) {
        if (value <= 0L) {
            stop(paste("Parameter", param$name, "must be positive"))
        }
    }
    if (param$nonnegative) {
        if (value < 0L) {
            stop(paste("Parameter", param$name, "must be non-negative"))
        }
    }
    value
}

check_parameter.float_parameter <- function (param, value) {
    value <- as.numeric(value)
    if (is.na(value)) {
        stop(paste("Parameter", param$name, "mustn't be NA"))
    }
    if (param$positive) {
        if (value <= 0) {
            stop(paste("Parameter", param$name, "must be positive"))
        }
    }
    if (param$nonnegative) {
        if (value < 0) {
            stop(paste("Parameter", param$name, "must be non-negative"))
        }
    }
    value
}

check_parameter.logical_parameter <- function (param, value) {
    value <- as.logical(value)
    if (is.na(value)) {
        stop(paste("Parameter", param$name, "mustn't be NA"))
    }
    value
}

check_parameter.mc_parameter <- function (param, value) {
    match.arg(value, param$mc)
}

check_parameter.file_parameter <- function(param, value) {
    file <- param$file
    if (!file.exists(file)) {
        stop(paste("File", file, "does not exist"))
    }
    file
}

check_parameter.default <- function(param, value) {
    stop("Invalid state")
}
