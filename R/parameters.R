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
                   mc = mc,
                   file = file), class = c(type, parameter_class))
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
        stop(paste("Parameter", param$name, "mustn't be zero length vector"))
    }
    UseMethod("check_parameter")
}

check_parameter.integer <- function (param, value) {
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
}

check_parameter.float <- function (param, value) {
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
}

check_parameter.logical <- function (param, value) {
    value <- as.logical(value)
    if (is.na(value)) {
        stop(paste("Parameter", param$name, "mustn't be NA"))
    }
}

check_parameter.mc <- function (param, value) {
    match.arg(value, param$mc)
}

check_parameter.file <- function(param, value) {
    file <- param$file
    if (!file.exists(file)) {
        stop(paste("File", file, "does not exist"))
    }
}

check_parameter.defaule <- function(param, value) {
    stop("Invalid state")
}
