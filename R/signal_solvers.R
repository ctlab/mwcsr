normalize_signals <- function(instance, signals) {
    if (is.null(signals)) {
        instance <- check_attr(instance, "weight", 0)
        V(instance)$signal <- paste0("V", 1:igraph::vcount(instance))
        E(instance)$signal <- paste0("E", 1:igraph::ecount(instance))
        signals <- setNames(c(V(instance)$weight, E(instance)$weight),
                            c(V(instance)$signal, E(instance)$signal))
    }

    instance <- check_attr(instance, "signal", NA, hint = "a default signal with zero weight.")
    observed_signals <- setdiff(union(V(instance)$signal, E(instance)$signal), NA)
    non_assigned_signals <- setdiff(observed_signals, names(signals))
    if (length(non_assigned_signals) != 0) {
        warning(paste0("No weight assigned for the following signals: ",
                       paste0(non_assigned_signals, collapse = ", ")))
        V(instance)$signal[V(instance)$signal %in% non_assigned_signals] <- NA
        E(instance)$signal[E(instance)$signal %in% non_assigned_signals] <- NA
    }

    # Enumerating signals
    map <- stats::setNames(1:(length(signals) + 1),
                           c(names(signals), "NA"))
    mapSignal <- function(x) {
        sapply(x, function(x) ifelse(is.na(x), map[length(map)], map[x]))
    }
    V(instance)$signal <- mapSignal(V(instance)$signal)
    E(instance)$signal <- mapSignal(E(instance)$signal)
    names(signals) <- mapSignal(names(signals))
    if (any(is.na(V(instance)$signal)) || any(is.na(E(instance)$signal))) {
        signals <- c(signals, setNames(0, map[length(map)]))
    }

    V(instance)$id <- 1:length(V(instance))
    E(instance)$id <- 1:length(E(instance))
    list(graph = instance, signals = signals)
}
