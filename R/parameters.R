parameter <- function(name,
                      type = c("integer", "float", "logical", "mc", "file", "char"),
                      is_null_possible = FALSE,
                      positive = FALSE,
                      nonnegative = FALSE,
                      mc = NULL) {
    structure(list(name = name,
                   type = type,
                   is_null_possible = is_null_possible,
                   positive = positive,
                   nonnegative = nonnegative,
                   mc = mc),
              class = c(paste0(type, "_parameter"), parameter_class))
}

params <- function(...) {
    l <- list(...)
    names(l) <- lapply(l, function(x) x$name)
    structure(l, class = "mwcs_solver_params")
}

#' @export
print.mwcs_solver_params <- function (x, ...) {
    df <- data.frame(name = sapply(x, function(x) x$name))
    df[["type"]] <- sapply(x, function(x) x$type)
    df[["value"]] <- sapply(x, function(x) {
        if (x$type == "mc") {
            quote <- function(x) paste0("\"", x, "\"")
            paste0("{", paste0(sapply(x$mc, quote), collapse = ", "), "}")
        } else if (x$positive) {
            "positive"
        } else if (x$nonnegative) {
            "non-negative"
        } else {
            ""
        }
    })
    print(df, row.names = FALSE, right = FALSE)
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

check_parameter.char_parameter <- function(param, value) {
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
    file <- normalizePath(value)
    if (!file.exists(file)) {
        stop(paste("File", file, "does not exist"))
    }
    file
}

check_parameter.default <- function(param, value) {
    stop("Invalid state")
}
