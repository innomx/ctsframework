#********************************************************************************
#
#     CTSFRAMEWORK
#
#********************************************************************************

cts_eval <- function(module, context, scenario, .scenario.path) {
    if (!is.null(seed <- attr(module, "seed"))) {
        save.seed <- .Random.seed
        on.exit({.Random.seed <- save.seed})
        set.seed(seed)
    }

    context2 <- list()
    context2$get_name <- function(x=module) get_name(x)
    context2$get_desc <- function(x=module) get_desc(x)

    if (!is.null(module)) {

        for (i in 1:length(module)) {
            nm <- names(module)[i]
            if (is.null(nm) || nm == "") {
                nm <- attr(module, "module")
            }
            context[[nm]] <- rlang::eval_tidy(module[[i]], c(context, context2))
        }
    }

    if (!is.null(outputs <- attr(module, "outputs"))) {
        .module.path <- file.path(.scenario.path, attr(module, "module"))
        dir.create(.module.path, showWarnings=FALSE, recursive=TRUE)
        outcontext <- c(context, list(.module.path=.module.path))
        outres <- rlang::eval_tidy(outputs, c(outcontext, context2))
        if (!is.null(outres)) {
            cts_save(outres, .path=.module.path)
        }
    }
    context
}

do_sim <- function(scenario, timestamp=format(Sys.time(), "%Y-%m-%d"), .path=file.path("cts_results", timestamp), .save=TRUE) {
    if (!is.null(seed <- attr(scenario, "seed"))) {
        save.seed <- .Random.seed
        on.exit({.Random.seed <- save.seed})
        set.seed(seed)
    }

    .scenario.path <- file.path(.path, attr(scenario, "name"))
    if (.save) {
        dir.create(.scenario.path, showWarnings=FALSE, recursive=TRUE)
        #saveRDS(scenario, file=file.path(.scenario.path, "scenario.rds"))
        cat("", file=file.path(.scenario.path, "scenario.R"))
        con <- file(file.path(.scenario.path, "scenario.R"), open="a")
        on.exit({close(con)})
        write_code(scenario, con=con)
    }

    file <- file.path(.path, "scenarios.csv")
    append <- file.exists(file) && file.size(file) > 0
    write.table(file=file, col.names=!append, row.names=FALSE, na="", sep=",", append=append,
        data.frame(
            timestamp     = Sys.time(),
            scenario      = get_name(scenario)               %||% NA,
            seed          = get_seed(scenario)               %||% NA,
            parameters    = get_name(scenario$parameters)    %||% NA,
            population    = get_name(scenario$population)    %||% NA,
            subject_model = get_name(scenario$subject_model) %||% NA,
            treatment     = get_name(scenario$treatment)     %||% NA,
            pkpd_model    = get_name(scenario$pkpd_model)    %||% NA,
            simulation    = get_name(scenario$simulation)    %||% NA,
            endpoints     = get_name(scenario$endpoints)     %||% NA,
            summarization = get_name(scenario$summarization) %||% NA
        )
    )
 
    x <- scenario

    f <- function(y) cts_eval(y, x, scenario, .scenario.path=.scenario.path)

    x <- f(scenario$parameters)
    x <- f(scenario$population)
    x <- f(scenario$subject_model)
    x <- f(scenario$treatment)
    x <- f(scenario$pkpd_model)
    x <- f(scenario$simulation)
    x <- f(scenario$endpoints)
    x <- f(scenario$summarization)

    x$scenario <- scenario

    invisible(structure(x, class="cts_scenario_evaluated", timestamp=Sys.time()))
}


#' Perform a clinical trial simulation
#'
#' @param x An object (e.g. a [scenario] or [scenario_list]).
#' @param ... Further arguments, ignored.
#' @export
CTS <- function(x, ...) UseMethod("CTS")

#' @export
CTS.cts_scenario <- function(x, timestamp=format(Sys.time(), "%Y-%m-%d"), .path=file.path("cts_results", timestamp), .save=FALSE) {
    res <- do_sim(x, timestamp=timestamp, .path=.path, .save=.save)
    fixup_scenario_table(.path=.path)
    invisible(res)
}

#' @export
CTS.cts_scenario_list <- function(x, timestamp=format(Sys.time(), "%Y-%m-%d"), .path=file.path("cts_results", timestamp), .save=FALSE) {
    res <- lapply(x, do_sim, timestamp=timestamp, .path=.path, .save=.save)
    fixup_scenario_table(.path=.path)
    names(res) <- sapply(res, get_name)
    class(res) <- "cts_scenario_evaluated_list" 
    invisible(res)
}

fixup_scenario_table <- function(.path, scenario_table="scenarios.csv") {
    x <- fread(file.path(.path, scenario_table), fill=T)
    x <- x[rev(!duplicated(rev(scenario)))]
    x <- x[scenario %in% dir(.path)]
    keep <- sapply(x, function(xx) !all(is.na(xx)))
    x <- x[, keep, with=F]
    fwrite(x, file.path(.path, scenario_table))
}


#' @export
replication <- function(scenario, nreplicates, seed=NULL, collate=NULL) {

    if (!is.null(attr(scenario, "seed"))) {
        warning("The scenario being replicated has a random seed. This is generally not desirable because it usually implies that the replicates will all be the same. Set the seed on the replication object instead.")
    }

    #do.call(scenario_list, c(replicate(nreplicates, scenario, simplify=FALSE), list(seed=seed, collate=collate)))
    as_scenario_list(replicate(nreplicates, scenario, simplify=FALSE), seed=seed, collate=collate)
}

#' @export
parameter_uncertainty <- function(scenario, ...) {
    stop("to be implemented")
}


#' @export
collate <- function(x, path, ..., outputs=NULL) UseMethod("collate")


#' @export
collate.cts_scenario_evaluated_list <- function(x, path, scenario_table=get_scenario_table(x), ..., outputs=NULL) {
    dat <- data.table(
        scenario = sapply(x, get_name)
    )[, getElement(x[[scenario]], name=path), by=scenario]
    dat <- setDT(scenario_table)[dat, on="scenario"]

    if (!is.null(outputs)) {
        outres <- rlang::eval_tidy(outputs, list(data=dat))
        if (!is.null(outres)) {
            cts_save(outres, .path=x)
        }
    }
    dat
}


#' @import data.table
#' @export
collate.character <- function(x, path, scenario_table="scenarios.csv", ..., outputs=NULL) {

    scenario_table <- fread(file.path(x, scenario_table))
    scenario_table <- scenario_table[rev(!duplicated(rev(scenario)))]

    scenario_files <- scenario_table[, .(file=file.path(x, scenario, path)), by=scenario]

    f <- function(file, scenario) {
        dat <- fread(file)
        data.table(scenario, dat)
    }
    dat <- mapply(f, scenario_files$file, scenario_files$scenario, SIMPLIFY=F)
    dat <- rbindlist(dat, fill=T)

    dat <- scenario_table[dat, on="scenario"]

    if (!is.null(outputs)) {
        outres <- rlang::eval_tidy(outputs, list(data=dat))
        if (!is.null(outres)) {
            cts_save(outres, .path=x)
        }
    }
    dat
}


#' @export
get_name <- function(x) {
    x %@% name
}


#' @export
get_desc <- function(x) {
    x %@% desc
}


#' @export
get_seed <- function(x) {
    x %@% seed
}


#' @export
get_scenario_table <- function(scenario_list) {
    x <- do.call(rbind, lapply(scenario_list, function(scenario) {
        data.frame(
            scenario      = get_name(scenario               ) %||% NA,
            seed          = get_seed(scenario               ) %||% NA,
            parameters    = get_name(scenario$parameters    ) %||% NA,
            population    = get_name(scenario$population    ) %||% NA,
            subject_model = get_name(scenario$subject_model ) %||% NA,
            treatment     = get_name(scenario$treatment     ) %||% NA,
            pkpd_model    = get_name(scenario$pkpd_model    ) %||% NA,
            simulation    = get_name(scenario$simulation    ) %||% NA,
            endpoints     = get_name(scenario$endpoints     ) %||% NA,
            summarization = get_name(scenario$summarization ) %||% NA
        )
    }))

    keep <- sapply(x, function(xx) !all(is.na(xx)))
    x[, keep, drop=F]
}


#********************************************************************************
#
#     MODULES
#
#********************************************************************************


#' CTS modules
#'
#' Modules are the components that make up a CTS scenario. The power of the
#' framework comes from the ability to mix and match modules to easily create a
#' multitute of scenarios.
#'
#' Each module has the same structure: it consists of a named list of
#' expressions that are run in sequence when the scenario is executed.
#'
#' @param ... A named list of expressions.
#' @param outputs (optional) A [cts_output] or [cts_output_list].
#' @param name Modules can be given a name.
#' @param seed Modules can be given a RNG seed, for reproducibility.
#' @return A module object.
#' @name cts_module
NULL


# Internal
make_cts_module <- function(module) {
    function(..., outputs=NULL, name=NULL, desc=NULL, seed=NULL, collate=NULL) {
        code    <- rlang::enquos(...)
        outputs <- rlang::enquo(outputs)
        if (rlang::quo_is_null(outputs)) outputs <- NULL
        counter <- next_counter(module)
        #if (is.null(name)) name <- paste0(module, counter)
        #if (is.null(name)) name <- ""
        structure(code, class=c(paste0("cts_", module), "cts_module"), module=module, counter=counter, outputs=outputs, name=name, desc=desc, seed=seed, collate=collate)
    }
}



#' @rdname cts_module
#' @export
parameters <- make_cts_module("parameters")

#' @rdname cts_module
#' @export
population <- make_cts_module("population")

#' @rdname cts_module
#' @export
subject_model  <- make_cts_module("subject_model")

#' @rdname cts_module
#' @export
treatment  <- make_cts_module("treatment")

#' @rdname cts_module
#' @export
pkpd_model <- make_cts_module("pkpd_model")

#' @rdname cts_module
#' @export
simulation <- make_cts_module("simulation")

#' @rdname cts_module
#' @export
endpoints <- make_cts_module("endpoints")

#' @rdname cts_module
#' @export
summarization <- make_cts_module("summarization")

#' CTS module (pipe interface)
#' @name cts_module_
NULL

#' @rdname cts_module_
#' @export
parameters_ <- function(scenario, ...) update(scenario, parameters=parameters(...))

#' @rdname cts_module_
#' @export
population_ <- function(scenario, ...) update(scenario, population=population(...))

#' @rdname cts_module_
#' @export
subject_model_  <- function(scenario, ...) update(scenario, subject_model=subject_model(...))

#' @rdname cts_module_
#' @export
treatment_  <- function(scenario, ...) update(scenario, treatment=treatment(...))

#' @rdname cts_module_
#' @export
pkpd_model_ <- function(scenario, ...) update(scenario, pkpd_model=pkpd_model(...))

#' @rdname cts_module_
#' @export
simulation_ <- function(scenario, ...) update(scenario, simulation=simulation(...))

#' @rdname cts_module_
#' @export
endpoints_ <- function(scenario, ...) update(scenario, endpoints=endpoints(...))

#' @rdname cts_module_
#' @export
summarization_ <- function(scenario, ...) update(scenario, summarization=summarization(...))




make_cts_module_grid <- function(module) {
    function(...) {
        f <- function(y) {
            #c(0, seq_along(y))
            c(seq_along(y))
        }

        args <- list(...)
        grid <- do.call(expand.grid, c(lapply(args, f), list(KEEP.OUT.ATTRS=FALSE)))

        lapply(split(grid, 1:nrow(grid)), function(grid.i) {
            a <- lapply(seq_along(grid.i)[grid.i > 0], function(j) args[[j]][[grid.i[[j]]]])
            names(a) <- names(args)[grid.i > 0]
            do.call(module, a)
        })
    }
}

#' CTS module grid
#' @name cts_module_grid
NULL

#' @rdname cts_module_grid
#' @export
parameters_grid <- make_cts_module_grid("parameters")

#' @rdname cts_module_grid
#' @export
population_grid <- make_cts_module_grid("population")

#' @rdname cts_module_grid
#' @export
subject_model_grid  <- make_cts_module_grid("subject_model")

#' @rdname cts_module_grid
#' @export
treatment_grid  <- make_cts_module_grid("treatment")

#' @rdname cts_module_grid
#' @export
pkpd_model_grid <- make_cts_module_grid("pkpd_model")

#' @rdname cts_module_grid
#' @export
simulation_grid <- make_cts_module_grid("simulation")

#' @rdname cts_module_grid
#' @export
endpoints_grid <- make_cts_module_grid("endpoints")

#' @rdname cts_module_grid
#' @export
summarization_grid <- make_cts_module_grid("summarization")



#********************************************************************************
#
#     SCENARIO
#
#********************************************************************************


#' Define a CTS scenario
#'
#' The scenario is the central concept of the CTS framework. A scenario is
#' defined by one or more modules. It is the ability to mix and match modules
#' to create a multitute of scenarios that gives the framework its power.
#'
#' @param parameters [parameters] module (optional).
#' @param population [population] module (optional).
#' @param subject_model [subject_model] module (optional).
#' @param treatment [treatment] module (optional).
#' @param pkpd_model [pkpd_model] module (optional).
#' @param simulation [simulation] module (optional).
#' @param endpoints [endpoints] module (optional).
#' @param summarization [summarization] module (optional).
#' @param ... Additional arguments, ignored.
#' @param name Optional scenario name.
#' @param desc Optional scenario description.
#' @param seed RNG seed, for reproducibility.
#' @param .silent Set to `TRUE` to supress informative messages.
#' @export
scenario <- function(
    parameters    = NULL,
    population    = NULL,
    subject_model = NULL,
    treatment     = NULL,
    pkpd_model    = NULL,
    simulation    = NULL,
    endpoints     = NULL,
    summarization = NULL,
    ...,
    name=NULL,
    desc=NULL,
    seed=make_seed(),
    .silent=FALSE
) {

    if (is.null(name)) {
        name <- paste0("Scenario", next_counter("cts.scenario.counter"))
    }

    if (!isTRUE(.silent)) {
        message(sprintf("Creating a new scenario named %s", name))
    }

    structure(
        list(
            parameters    = parameters,
            population    = population,
            subject_model = subject_model,
            treatment     = treatment,
            pkpd_model    = pkpd_model,
            simulation    = simulation,
            endpoints     = endpoints,
            summarization = summarization
        ), class="cts_scenario", name=name, desc=desc, seed=seed)
}



#' Derive a new scenario from an existing one
#'
#' A new scenario can be derived from an existing one by overriding any of its
#' modules, or specific components of its modules.
#'
#' @param base_scenario A [scenario] which serves as the base scenario from which to derive.
#' @param ... Further argument passed to [scenario()].
#' @inheritParams scenario
#' @export
derive_scenario <- function(

    base_scenario,

    parameters    = NULL,
    population    = NULL,
    subject_model = NULL,
    treatment     = NULL,
    pkpd_model    = NULL,
    simulation    = NULL,
    endpoints     = NULL,
    summarization = NULL,
    ...,
    name=NULL,
    desc=NULL,
    seed=NULL
) {

    #if (is.null(seed) && !is.null(attr(base_scenario, "seed"))) {
    #    warning("The derived scenario has the same random seed as the base scenario.")
    #}

    x <- scenario(
        parameters    = merge_lists(parameters    , base_scenario$parameters),
        population    = merge_lists(population    , base_scenario$population),
        subject_model = merge_lists(subject_model , base_scenario$subject_model),
        treatment     = merge_lists(treatment     , base_scenario$treatment),
        pkpd_model    = merge_lists(pkpd_model    , base_scenario$pkpd_model),
        simulation    = merge_lists(simulation    , base_scenario$simulation),
        endpoints     = merge_lists(endpoints     , base_scenario$endpoints),
        summarization = merge_lists(summarization , base_scenario$summarization),
        ...,
        name=name,
        desc=desc,
        seed=seed %||% attr(base_scenario, "seed")
    )

    for (a in c("name", "desc", "seed")) {
        attr(x$parameters    , a) <- attr(parameters    , a, exact=TRUE) %||% attr(base_scenario$parameters    , a, exact=TRUE)
        attr(x$population    , a) <- attr(population    , a, exact=TRUE) %||% attr(base_scenario$population    , a, exact=TRUE)
        attr(x$subject_model , a) <- attr(subject_model , a, exact=TRUE) %||% attr(base_scenario$subject_model , a, exact=TRUE)
        attr(x$treatment     , a) <- attr(treatment     , a, exact=TRUE) %||% attr(base_scenario$treatment     , a, exact=TRUE)
        attr(x$pkpd_model    , a) <- attr(pkpd_model    , a, exact=TRUE) %||% attr(base_scenario$pkpd_model    , a, exact=TRUE)
        attr(x$simulation    , a) <- attr(simulation    , a, exact=TRUE) %||% attr(base_scenario$simulation    , a, exact=TRUE)
        attr(x$endpoints     , a) <- attr(endpoints     , a, exact=TRUE) %||% attr(base_scenario$endpoints     , a, exact=TRUE)
        attr(x$summarization , a) <- attr(summarization , a, exact=TRUE) %||% attr(base_scenario$summarization , a, exact=TRUE)
    }

    x
}

#' Update an existing scenario
#'
#' @param object The [scenario] to be updated.
#' @param ... Further arguments passed to [derive_scenario()].
#' @inheritParams derive_scenario
#' @inheritParams scenario
#' @export
update.cts_scenario <- function(object, ..., .silent=TRUE) {
    derive_scenario(object, ...,
        name=get_name(object),
        #desc=attr(object, "desc"),
        .silent=.silent
    )
}

#' Create a scenario list
#'
#' Multiples scenarios are combined into a list that can be run as a unit.
#'
#' @param ... Any number of [scenario]s.
#' @export
scenario_list <- function(..., collate=NULL) {
    collate <- rlang::enquo(collate)
    if (rlang::quo_is_null(collate)) collate <- NULL
    if (!all(sapply(list(...), inherits, what="cts_scenario"))) {
        stop("A scenario_list should only contain scenarios.")
    }
    structure(list(...), class=c("cts_scenario_list"), collate=collate)
}

as_scenario_list <- function(x, ...) {
    do.call(scenario_list, c(x, list(...)))
    #structure(x, class=c("cts_scenario_list"), ...)
}

#' Create a scenario grid
#'
#' A scenario grid is like a [scenario_list] except that all combinations of
#' multiple different factors are generated.
#'
#' @inheritParams derive_scenario
#' @inheritParams scenario
#' @export
scenario_grid <- function(

    base_scenario,

    parameters    = NULL,
    population    = NULL,
    subject_model = NULL,
    treatment     = NULL,
    pkpd_model    = NULL,
    simulation    = NULL,
    endpoints     = NULL,
    summarization = NULL,
    ...
) {

    f <- function(y) {
        #c(0, seq_along(y))
        c(seq_along(y))
    }

    args <- list(
        parameters    = parameters,
        population    = population,
        subject_model = subject_model,
        treatment     = treatment,
        pkpd_model    = pkpd_model,
        simulation    = simulation,
        endpoints     = endpoints,
        summarization = summarization
    )
    args <- args[!sapply(args, is.null)]

    grid <- do.call(expand.grid, c(lapply(args, f), list(KEEP.OUT.ATTRS=FALSE)))

    sclist <- lapply(split(grid, 1:nrow(grid)), function(grid.i) {
        a <- lapply(seq_along(grid.i)[grid.i > 0], function(j) args[[j]][[grid.i[[j]]]])
        names(a) <- names(args)[grid.i > 0]
        do.call(derive_scenario, c(list(base_scenario), a))
    })

    #do.call(scenario_list, c(sclist, list(...)))
    as_scenario_list(sclist, ...)
}

# print.cts_scenario_grid <- function(x, ...) {
#     print(as.data.frame(x)[, apply(x, 2, sum) > 0, drop=F])
# }



#********************************************************************************
#
#     WRITE_CODE
#
#********************************************************************************


#' @export
write_code <- function(x, ..., con) UseMethod("write_code")

#' @export
write_code.NULL <- function(x, ..., con) {
    invisible(NULL) 
}

#' @export
write_code.cts_scenario <- function(x, ..., con) {
    cat("save.seed <- list()\n\n", file=con)
    if (!is.null(seed <- attr(x, "seed"))) {
        cat("save.seed <- push(save.seed, .Random.seed)\n", file=con)
        cat("set.seed(", seed, ")\n", file=con)
    }

    lapply(x, write_code, ..., con=con)

    if (!is.null(seed <- attr(x, "seed"))) {
        cat("\n.Random.seed <- top(save.seed)\n", file=con)
        cat("save.seed <- pop(save.seed)\n", file=con)
    }
    invisible(NULL) 
}

#' @export
write_code.cts_module <- function(x, ..., con) {
    module.name <- attr(x, "module")
    cat("\n\n### ", module.name, "\n\n", file=con)

    if (!is.null(seed <- attr(x, "seed"))) {
        cat("save.seed <- push(save.seed, .Random.seed)\n", file=con)
        cat("set.seed(", seed, ")\n", file=con)
    }

    for (i in 1:length(x)) {
        if (!is.null(names(x))) nm <- names(x)[i]
        if (is.null(names(x)) || nm=="") nm <- module.name
        cat(nm, " <- ", file=con)
        cat(deparse1(rlang::quo_get_expr(x[[i]]), width.cutoff=60L, collapse="\n"), "\n\n", file=con)
    }

    # Outputs will not be generated
    #if (!is.null(outputs <- attr(x, "outputs"))) {
    #    cat("outputs <- ", file=con)
    #    cat(deparse1(rlang::quo_get_expr(outputs), width.cutoff=60L, collapse="\n"), "\n\n", file=con)
    #}

    if (!is.null(seed <- attr(x, "seed"))) {
        cat("\n.Random.seed <- top(save.seed)\n", file=con)
        cat("save.seed <- pop(save.seed)\n", file=con)
    }
    invisible(NULL) 
}




#********************************************************************************
#
#     OUTPUTS
#
#********************************************************************************


#' @export
output_list <- function(...) {
    structure(list(...), class="cts_output_list")
}

#' @export
graphical_output <- function(x, ..., name, width, height) {
    structure(x, class=c("cts_graphical_output", "cts_output", class(x)), name=name, width=width, height=height)
}

#' @export
tabular_output <- function(x, ..., name) {
    structure(x, class=c("cts_tabular_output", "cts_output", class(x)), name=name)
}

#' @export
csv_output <- function(x, ..., name) {
    structure(x, class=c("cts_csv_output", "cts_output", class(x)), name=name)
}

#' @export
yaml_output <- function(x, ..., name) {
    structure(x, class=c("cts_yaml_output", "cts_output", class(x)), name=name)
}

#' @export
rds_output <- function(x, ..., name) {
    structure(x, class=c("cts_rds_output", "cts_output", class(x)), name=name)
}

#' @export
cts_save <- function(x, ...) UseMethod("cts_save")

#' @export
cts_save.cts_output_list <- function(x, ..., .path) {
    if (is.null(names(x))) {
        names(x) <- paste0("output", 1:length(x))
    }
    for (i in 1:length(x)) {
        x[[i]] %@% name <- get_name(x[[i]]) %||% names(x)[i]
    }
    lapply(x, cts_save, .path=.path)
}

#' @export
cts_save.cts_graphical_output <- function(x, ..., .path) {
    ggplot2::ggsave(plot=x, filename=file.path(.path, paste0(get_name(x), ".png")),
        width=(x %@% width), height=(x %@% height), ...)
}

#' @export
cts_save.cts_tabular_output <- function(x, ..., .path) {
}

#' @export
cts_save.cts_csv_output <- function(x, ..., .path) {
    data.table::fwrite(x=x, file=file.path(.path, paste0(get_name(x), ".csv")), ...)
}

#' @export
cts_save.cts_yaml_output <- function(x, ..., .path) {
}

#' @export
cts_save.cts_rds_output <- function(x, ..., .path) {
}



#********************************************************************************
#
#     UTILS
#
#********************************************************************************



# Functions copied from rlang
`%||%`  <- rlang::`%||%`
`%@%`   <- rlang::`%@%`
`%@%<-` <- rlang::`%@%<-`

merge_lists <- function(a, b) { b[names(a)] <- a ; b }

make_seed <- function() {
    round(1e8 * runif(1))
}

closure_counter <- (function() {
    myenv <- new.env(parent=emptyenv())
    list(
        next_counter = function(counter) {
            stopifnot(is.character(counter) && length(counter) == 1)
            n <- get0(counter, envir=myenv, mode="integer", ifnotfound=1)
            assign(counter, (n+1), envir=myenv)
            n
        },
        set_counter = function(counter, value) {
            stopifnot(is.character(counter) && length(counter) == 1)
            assign(counter, value,envir=myenv)
        },
        reset_counter = function(counter, .all=FALSE) {
            if (.all) {
                rm(list=ls(myenv), envir=myenv)
            } else {
                stopifnot(is.character(counter))
                rm(list=intersect(counter, ls(myenv)),envir=myenv)
            }
        },
        list_counters = function() {
            sapply(ls(myenv), get, envir=myenv)
        }
    )
})()

next_counter  <- closure_counter$next_counter
set_counter   <- closure_counter$set_counter
reset_counter <- closure_counter$reset_counter
list_counters <- closure_counter$list_counters


#' Quantile interval
#'
#' A level $α$ quantile interval of `x` is the interval whose endpoints are the
#' upper and lower $(1-α)/2$ quantiles of `x`. In other words, it covers a
#' $(1-α)$ fraction of the data, with an equal number of points outside the
#' interval in each tail.
#'
#' @param x A numeric vector.
#' @param level The level of the interval, typically 0.95, 0.9 or 0.5 (which
#' corresponds to the interquartile range).  Determines the quantiles that
#' correspond to the endpoints of the interval.  Can also be a vector to
#' compute multiple intervals at once.
#' @param type The quantile algorithm (see [stats::quantile()]).
#' @return A named list of numbers, giving the lower and upper limits of the
#' interval, as well as the median.
#' @examples
#' x <- 1:1000
#' qi(x, 0.9)
#' qi(x, c(0.95, 0.9, 0.5))
#' @export
qi <- function(x, level=0.95, type=7) {
    # XXX check inputs
    probs <- sort(unique(c((1 - level)/2, 0.5, 1 - (1 - level)/2)))
    qname <- paste0("q", probs)
    q <- quantile(x, probs=probs, type=type, names=FALSE, na.rm=TRUE)
    setNames(as.list(q), qname)
}


#********************************************************************************
#
#     STACK-FUNCTIONS
#
#********************************************************************************


#' Stack Functions
#'
#' A simple stack implementation. Provides the basic functions `push`, `pop`
#' and `top`. The stack itself is implemented as a simple list.
#'
#' @details In the event the stack is empty, both `pop` and `top` issue a warning.
#'
#' @param stack A list.
#' @param value Any value or object.
#' @return `push` and `pop` return the stack in its state after the specified
#' operation. `top` return the value currently at the top of the stack, or null
#' if the stack is empty.
#' @examples
#' stack <- list()
#' stack <- stack |> push(3) |> push(5) |> push(7)
#' top(stack)
#' stack <- pop(stack)
#' top(stack)
#' stack <- pop(stack)
#' top(stack)
#' stack <- pop(stack)
#' top(stack)
#' stack <- pop(stack)
#' @name stack-functions
NULL

#' @rdname stack-functions
#' @export
push <- function(stack, value) {
    stack <- c(list(value), stack)
    invisible(stack)
}

#' @rdname stack-functions
#' @export
pop <- function(stack) {
    if (length(stack) == 0) {
        warning("Stack is empty")
    } else {
        stack <- stack[-1]
    }
    invisible(stack)
}

#' @rdname stack-functions
#' @export
top <- function(stack) {
    if (length(stack) == 0) {
        warning("Stack is empty")
        return(NULL)
    }
    stack[[1]]
}

