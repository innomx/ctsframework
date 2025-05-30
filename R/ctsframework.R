#********************************************************************************
#
#     CTSFRAMEWORK
#
#********************************************************************************

# @internal
cts_eval <- function(module, context, scenario, .scenario.path) {
    if (!is.null(seed <- get_seed(module))) {
        save.seed <- .Random.seed
        on.exit({.Random.seed <- save.seed})
        set.seed(seed)
    }

    if (!is.null(module)) {

        for (i in seq_along(module)) {
            nm <- names(module)[i]
            if (is.null(nm) || nm == "") {
                nm <- get_module(module)
            }
            context[[nm]] <- rlang::eval_tidy(module[[i]], context)
        }
    }

    if (!is.null(outputs <- get_outputs(module))) {
        .module.path <- file.path(.scenario.path, get_module(module))
        dir.create(.module.path, showWarnings=FALSE, recursive=TRUE)
        outcontext <- c(context, list(.module.path=.module.path))
        outres <- rlang::eval_tidy(outputs, outcontext)
        if (!is.null(outres)) {
            cts_save(outres, .path=.module.path)
        }
    }
    context
}

# @internal
do_sim <- function(scenario, timestamp=Sys.time(), .path=file.path("cts_results", as.Date(timestamp)), .save=TRUE) {

    if (!is.null(seed <- get_seed(scenario))) {
        save.seed <- .Random.seed
        on.exit({.Random.seed <- save.seed})
        set.seed(seed)
    }

    .scenario.path <- file.path(.path, get_name(scenario))

    if (.save) {
        dir.create(.scenario.path, showWarnings=FALSE, recursive=TRUE)
        #saveRDS(scenario, file=file.path(.scenario.path, "scenario.rds"))
        cat("", file=file.path(.scenario.path, "scenario.R"))
        con <- file(file.path(.scenario.path, "scenario.R"), open="a")
        on.exit({close(con)})
        write_code(scenario, con=con)
    }

    x <- scenario

    f <- function(y) cts_eval(y, x, scenario, .scenario.path=.scenario.path)

    #x <- f(scenario$parameters)
    #x <- f(scenario$population)
    #x <- f(scenario$subject_model)
    #x <- f(scenario$treatment)
    #x <- f(scenario$pkpd_model)
    #x <- f(scenario$simulation)
    #x <- f(scenario$endpoints)
    #x <- f(scenario$summarization)

    for (m in names(scenario)) {
        x <- f(scenario[[m]])
        if (!is.null(scenario[[m]]$name)) {
            set_name(scenario[[m]]) <- x$name
        }
    }

    x$scenario <- scenario

    any_outputs <- !all(sapply(lapply(scenario, get_outputs), is.null))

    if (any_outputs) {
        dir.create(.path, showWarnings=FALSE, recursive=TRUE)
        file <- file.path(.path, "scenarios.csv")
        append <- file.exists(file) && file.size(file) > 0
        write.table(file=file, col.names=!append, row.names=FALSE, na="", sep=",", append=append,
            data.table(
                timestamp     = timestamp,
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
    }
 

    invisible(structure(x, class="cts_scenario_evaluated", timestamp=timestamp))
}


#' Perform a clinical trial simulation
#'
#' This is the function that will actually run the [scenario] or
#' [scenario_list].
#'
#' @param x An object (e.g. a [scenario] or [scenario_list]).
#' @param ... Additional arguments (currently ignored).
#' @param timestamp This gets used to record when the simulation was run. It is
#' also used (by default) in the construction of the path containing the
#' simulation results (the date portion only). This helps to prevent overwriting
#' previous results by accident.
#' @param .path Path where the simulation results will be saved. Typically, a
#' folder structure will be created under this directory, with a subdirectory
#' for each scenario. Within each scenario directory, outputs are saved under a
#' subdirectory named for the module in which it is created.  Several other
#' files will be created as well, such as `scenarios.csv` (a record of the
#' scenarios that have been run in this directory) and a `scenario.R` for each
#' scenario (if `.save = TRUE`).
#' @param .save Logical. Should the `scenario.R` file containing the code
#' that can basically be run to re-create the simulated scenario be generated?
#' The point is not really to run this code, but it is nice to have as a record
#' of what was done to generate the outputs.
#' @return An object of class `cts_scenario_evaluated` or
#' `cts_scenario_evaluated_list` depending on whether a single [scenario] or a
#' [scenario_list] was passed in.
#' @examples
#' # A silly example that does nothing. For realistic examples, see the vignette.
#' CTS(scenario(name="foo"))
#' @export
CTS <- function(x, ...) UseMethod("CTS")

#' @rdname CTS
#' @export
CTS.cts_scenario <- function(x, ..., timestamp=Sys.time(), .path=file.path("cts_results", as.Date(timestamp)), .save=FALSE) {
    res <- do_sim(x, timestamp=timestamp, .path=.path, .save=.save)
    fixup_scenario_table(.path=.path)
    invisible(res)
}

#' @rdname CTS
#' @export
CTS.cts_scenario_list <- function(x, ..., timestamp=Sys.time(), .path=file.path("cts_results", as.Date(timestamp)), .save=FALSE) {
    res <- lapply(x, do_sim, timestamp=timestamp, .path=.path, .save=.save)
    fixup_scenario_table(.path=.path)
    names(res) <- sapply(res, get_name)
    class(res) <- "cts_scenario_evaluated_list" 
    invisible(res)
}

# @internal
fixup_scenario_table <- function(.path, scenario_table="scenarios.csv") {
    if (file.exists(file.path(.path, scenario_table))) {
        x <- fread(file.path(.path, scenario_table), fill=T)
        x <- x[rev(!duplicated(rev(scenario)))]
        x <- x[scenario %in% dir(.path)]
        if (nrow(x)==0 || ncol(x)==0) {
            file.remove(file.path(.path, scenario_table))
        } else {
            suppressWarnings(
                data.table::fwrite(x, file.path(.path, scenario_table))
            )
        }
    }
}


# #' Scenario replication
# #'  
# #' @param scenario The [scenario] to replicate.
# #' @param nreplicates The number of times to replicate the scenario.
# #' @export
# replication <- function(scenario, nreplicates) {
# 
#     stop("Not implemented yet")
# 
# }

# #' Add parameter uncertainty to a scenario
# #'
# #' @param scenario The [scenario] to add parameter uncertainty to.
# #' @param ... Additional arguments (currently ignored).
# #' @export
# parameter_uncertainty <- function(scenario, ...) {
# 
#     stop("Not implemented yet")
# 
# }


#' Collate simulation results from multiple scenarios
#'
#' This function can be used after running multiple scenarios with `CTS` to
#' combine the results for further processing and output generation.
#'
#' @param x An object.
#' @param ... Additional arguments (currently ignored).
#' @param outputs (optional) A [cts_output] or [cts_output_list].
#' @param .path The path containing the simulation results. This directory
#' should contain the file `scenarios.csv` as well as subfolders for each
#' scenario, as generated by [CTS()].
#' @return A `data.frame`.
#' @export
collate <- function(x, ..., outputs=NULL) UseMethod("collate")


# #' @export
# collate.cts_scenario_evaluated_list <- function(x, ..., what, scenario_table=get_scenario_table(x), outputs=NULL) {
#     dat <- data.table(
#         scenario = sapply(x, get_name)
#     )[, getElement(x[[scenario]], name=what), by=scenario]
#     dat <- setDT(scenario_table)[dat, on="scenario"]
# 
#     if (!is.null(outputs)) {
#         outres <- rlang::eval_tidy(outputs, list(data=dat))
#         if (!is.null(outres)) {
#             cts_save(outres, .path=x)
#         }
#     }
#     dat
# }


#' @rdname collate
#' @export
collate.character <- function(x, ..., outputs=NULL, .path=file.path("cts_results", as.Date(Sys.time()))) {

    scenario_table <- data.table::fread(file.path(.path, "scenarios.csv"))
    scenario_table <- scenario_table[rev(!duplicated(rev(scenario)))]

    # Remove empty columns
    keep <- sapply(scenario_table, function(x) !all(is.na(x)))
    scenario_table <- scenario_table[, keep, with=F]

    scenario_files <- scenario_table[, list(file=file.path(.path, scenario, x)), by=scenario]

    f <- function(file, scenario) {
        dat <- data.table::fread(file)
        data.table(scenario, dat)
    }
    dat <- mapply(f, scenario_files$file, scenario_files$scenario, SIMPLIFY=F)
    dat <- data.table::rbindlist(dat, fill=T)

    dat <- scenario_table[dat, on="scenario"]

    if (!is.null(outputs)) {
        outres <- rlang::eval_tidy(outputs, list(data=dat))
        if (!is.null(outres)) {
            cts_save(outres, .path=.path)
        }
    }
    dat
}


# @internal
get_name <- function(x) {
    x %@% "name"
}

# @internal
`set_name<-` <- function(x, value) {
    x %@% "name" <- value
    x
}

# @internal
get_seed <- function(x) {
    x %@% "seed"
}

# @internal
get_module <- function(x) {
    x %@% "module"
}

# @internal
get_outputs <- function(x) {
    x %@% "outputs"
}



# #' Get the scenario table for a list of scenarios
# #'
# #' @param scenario_list The [scenario_list] for which to get the table.
# #' @return A `data.frame`.
# #' @export
# get_scenario_table <- function(scenario_list) {
#     data.table::rbindlist(fill=T, lapply(scenario_list, function(scenario) {
#         data.table(
#             scenario      = get_name(scenario               ),
#             seed          = get_seed(scenario               ),
#             parameters    = get_name(scenario$parameters    ),
#             population    = get_name(scenario$population    ),
#             subject_model = get_name(scenario$subject_model ),
#             treatment     = get_name(scenario$treatment     ),
#             pkpd_model    = get_name(scenario$pkpd_model    ),
#             simulation    = get_name(scenario$simulation    ),
#             endpoints     = get_name(scenario$endpoints     ),
#             summarization = get_name(scenario$summarization )
#         )
#     })) |> as.data.frame()
# }


#********************************************************************************
#
#     MODULES
#
#********************************************************************************


#' CTS modules
#'
#' Modules are the components that make up a CTS scenario. The power of the
#' framework comes from the ability to mix and match modules to easily create a
#' multitude of scenarios.
#'
#' Each module has the same structure: it consists of a named list of
#' expressions that are run in sequence when the scenario is executed. These
#' functions that create modules are not meant to be called on their own, but
#' rather inside a [scenario].
#'
#' @param ... A named list of expressions.
#' @param outputs (optional) A [cts_output] or [cts_output_list].
#' @param seed (optional) Modules can be given a RNG seed, for reproducibility.
#' @return A `cts_module` object.
#' @examples
#' # A silly example that does nothing. For realistic examples, see the vignette.
#' CTS(scenario(name="foo",
#'     parameters    = parameters(),
#'     population    = population(),
#'     subject_model = subject_model(),
#'     treatment     = treatment(),
#'     pkpd_model    = pkpd_model(),
#'     simulation    = simulation(),
#'     endpoints     = endpoints(),
#'     summarization = summarization()
#' ))
#'
#' # There is also a syntax that uses pipes to build up the scenario
#' # (note the mandatory underscore '_' at the end of each module name):
#' scenario(name="foo") |>
#'     parameters_() |>
#'     population_() |>
#'     subject_model_() |>
#'     treatment_() |>
#'     pkpd_model_() |>
#'     simulation_() |>
#'     endpoints_() |>
#'     summarization_() |>
#'     CTS()
#'
#' @name cts_module
NULL


# @internal
make_cts_module <- function(module) {
    function(..., outputs=NULL, seed=NULL) {
        code    <- rlang::enquos(...)
        outputs <- rlang::enquo(outputs)
        if (rlang::quo_is_null(outputs)) outputs <- NULL
        counter <- next_counter(paste0("cts.", module))
        structure(code, class=c(paste0("cts_", module), "cts_module"), module=module, counter=counter, outputs=outputs, seed=seed)
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
#'
#' @param scenario A [scenario]. 
#' @param ... A named list of expressions.
#' @return A module object.
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




# @internal
make_cts_module_grid <- function(module) {
    function(...) {
        f <- function(y) {
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
#'
#' @param ... Any objects (typically vectors or lists).
#' @return A `list` of [cts_module] objects.
#' @examples
#' # This will create 5 scenarios, each with a different value for the parameter x:
#' scenario_grid(scenario(name="foo"), parameters=parameters_grid(x=1:5))
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


#' @exportS3Method base::print
print.cts_module <- function(x, ...) {
    cat(paste0(get_module(x), ": ", get_name(x), "\n"))
}


#********************************************************************************
#
#     SCENARIO
#
#********************************************************************************


#' Define a CTS scenario
#'
#' The scenario is the central concept of the CTS framework. A scenario is
#' defined by one or more modules. It is the ability to mix and match modules
#' to create a multitude of scenarios that gives the framework its power.
#'
#' @param parameters [parameters] module (optional).
#' @param population [population] module (optional).
#' @param subject_model [subject_model] module (optional).
#' @param treatment [treatment] module (optional).
#' @param pkpd_model [pkpd_model] module (optional).
#' @param simulation [simulation] module (optional).
#' @param endpoints [endpoints] module (optional).
#' @param summarization [summarization] module (optional).
#' @param ... Additional arguments (currently ignored).
#' @param name Optional scenario name. If omitted, a default name will be
#' generated such as `ScenarioX` where `X` is a number that is incremented for
#' each new scenario created.
#' @param seed RNG seed, for reproducibility.
#' @param .silent Set to `TRUE` to suppress informative messages.
#' @return A `cts_scenario` object.
#' @examples
#' # A silly example that does nothing. For realistic examples, see the vignette.
#' CTS(scenario(name="foo",
#'     parameters    = parameters(),
#'     population    = population(),
#'     subject_model = subject_model(),
#'     treatment     = treatment(),
#'     pkpd_model    = pkpd_model(),
#'     simulation    = simulation(),
#'     endpoints     = endpoints(),
#'     summarization = summarization()
#' ))
#'
#' # There is also a syntax that uses pipes to build up the scenario
#' # (note the mandatory underscore '_' at the end of each module name):
#' scenario(name="foo") |>
#'     parameters_() |>
#'     population_() |>
#'     subject_model_() |>
#'     treatment_() |>
#'     pkpd_model_() |>
#'     simulation_() |>
#'     endpoints_() |>
#'     summarization_() |>
#'     CTS()
#'
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
    name = NULL,
    seed = make_seed(),
    .silent = FALSE
) {

    if (is.null(name)) {
        name <- paste0("Scenario", next_counter("cts.scenario"))
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
        ), class="cts_scenario", name=name, seed=seed)
}



#' Derive a new scenario from an existing one
#'
#' A new scenario can be derived from an existing one by overriding any of its
#' modules, or specific components of its modules.
#'
#' @note By default, the derived scenario has the same RNG seed as the base
#' scenario. This can be changed by specifying the `seed` argument.
#'
#' @param base_scenario A [scenario] which serves as the base scenario from which to derive.
#' @param ... Further argument passed to [scenario()].
#' @inheritParams scenario
#' @return A `cts_scenario` object.
#' @examples
#' # A silly example that does nothing. For realistic examples, see the vignette.
#'
#' # Create a base scenario in which x=1
#' base <- scenario(name="base",
#'   parameters = parameters(x=1),
#'   simulation = simulation({cat("x =", x, "\n")})
#' )
#'
#' # Overwrite the value of x in the derived scenario
#' derived <- derive_scenario(base, name="derived", parameters=parameters(x=2))
#'
#' # Run both scenarios and observed the output
#' CTS(base)
#' CTS(derived)
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
    name = NULL,
    seed = NULL
) {

    x <- scenario(
        parameters    = merge_modules(parameters    , base_scenario$parameters),
        population    = merge_modules(population    , base_scenario$population),
        subject_model = merge_modules(subject_model , base_scenario$subject_model),
        treatment     = merge_modules(treatment     , base_scenario$treatment),
        pkpd_model    = merge_modules(pkpd_model    , base_scenario$pkpd_model),
        simulation    = merge_modules(simulation    , base_scenario$simulation),
        endpoints     = merge_modules(endpoints     , base_scenario$endpoints),
        summarization = merge_modules(summarization , base_scenario$summarization),
        ...,
        name=name,
        seed=seed %||% get_seed(base_scenario)
    )

    x
}

# Update an existing scenario
#
# @param object The [scenario] to be updated.
# @param ... Further arguments passed to [derive_scenario()].
# @inheritParams derive_scenario
# @inheritParams scenario
# @internal
#' @export
update.cts_scenario <- function(object, ..., .silent=TRUE) {
    derive_scenario(object, ...,
        name=get_name(object),
        .silent=.silent
    )
}

#' Create a scenario list
#'
#' Multiples scenarios are combined into a list that can be run as a unit.
#'
#' @param ... Any number of [scenario] objects.
#' @return A `cts_scenario_list` object.
#' @examples
#' # A silly example that does nothing. For realistic examples, see the vignette.
#' CTS(scenario_list(
#'     scenario(name="foo"),
#'     scenario(name="bar")
#' ))
#' @export
scenario_list <- function(...) {
    if (!all(sapply(list(...), inherits, what="cts_scenario"))) {
        stop("A scenario_list should only contain scenarios.")
    }
    sclist <- list(...)
    names(sclist) <- sapply(sclist, get_name)
    structure(sclist, class=c("cts_scenario_list"))
}

# @internal
as_scenario_list <- function(x, ...) {
    do.call(scenario_list, c(x, list(...)))
}

#' Create a scenario grid
#'
#' A scenario grid is like a [scenario_list] except that all combinations of
#' multiple different factors are generated.
#'
#' @inheritParams derive_scenario
#' @inheritParams scenario
#' @return A `cts_scenario_list` object.
#' @examples
#' # A silly example that does nothing. For realistic examples, see the vignette.
#'
#' # This will create 12 scenarios in total, representing all
#' # combinations of x, y and treatment.
#' CTS(scenario_grid(
#'     scenario(name="base", simulation=simulation({cat(x, y, treatment, "\n")})),
#'     parameters=parameters_grid(x=1:2, y=3:5),
#'     treatment=list(
#'         treatment("A"),
#'         treatment("B")
#'     )
#' ))
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
        do.call(derive_scenario, c(list(base_scenario=base_scenario), a))
    })

    names(sclist) <- sapply(sclist, get_name)

    as_scenario_list(sclist, ...)
}

#' @exportS3Method base::print
print.cts_scenario <- function(x, ...) {
    cat(paste0("Scenario: ", get_name(x), "\n"))
}

#' @exportS3Method base::print
print.cts_scenario_list <- function(x, ...) {
    cat(paste0("Scenario list:", "\n"))
    for (i in seq_along(x)) {
        cat(paste0("  ", get_name(x[[i]]), "\n"))
    }
}



#********************************************************************************
#
#     WRITE_CODE
#
#********************************************************************************


# @internal
write_code <- function(x, ..., con) UseMethod("write_code")

# @internal
#' @export
write_code.NULL <- function(x, ..., con) {
    invisible(NULL) 
}

# @internal
#' @export
write_code.cts_scenario <- function(x, ..., con) {
    cat("save.seed <- list()\n\n", file=con)
    if (!is.null(seed <- get_seed(x))) {
        cat("save.seed <- push(save.seed, .Random.seed)\n", file=con)
        cat("set.seed(", seed, ")\n", file=con)
    }

    lapply(x, write_code, ..., con=con)

    if (!is.null(seed <- get_seed(x))) {
        cat("\n.Random.seed <- top(save.seed)\n", file=con)
        cat("save.seed <- pop(save.seed)\n", file=con)
    }
    invisible(NULL) 
}

# @internal
#' @export
write_code.cts_module <- function(x, ..., con) {
    module.name <- get_module(x)
    cat("\n\n### ", module.name, "\n\n", file=con)

    if (!is.null(seed <- get_seed(x))) {
        cat("save.seed <- push(save.seed, .Random.seed)\n", file=con)
        cat("set.seed(", seed, ")\n", file=con)
    }

    for (i in seq_along(x)) {
        if (!is.null(names(x))) nm <- names(x)[i]
        if (is.null(names(x)) || nm=="") nm <- module.name
        cat(nm, " <- ", file=con)
        cat(deparse1(rlang::quo_get_expr(x[[i]]), width.cutoff=60L, collapse="\n"), "\n\n", file=con)
    }

    # Outputs will not be generated
    #if (!is.null(outputs <- get_outputs(x))) {
    #    cat("outputs <- ", file=con)
    #    cat(deparse1(rlang::quo_get_expr(outputs), width.cutoff=60L, collapse="\n"), "\n\n", file=con)
    #}

    if (!is.null(seed <- get_seed(x))) {
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

#' CTS outputs
#'
#' Modules are the components that make up a CTS scenario. The power of the
#' framework comes from the ability to mix and match modules to easily create a
#' multitude of scenarios.
#'
#' Each module has the same structure: it consists of a named list of
#' expressions that are run in sequence when the scenario is executed.
#'
#' @param x An expression.
#' @param name A name. It is used to construct the filename for the saved output.
#' @param width For graphical outputs, the width in inches.
#' @param height For graphical outputs, the height in inches.
#' @param ... Additional arguments (currently ignored).
#' @return An object that inherits from the class `cts_output`.
#' @name cts_output
NULL

#' CTS output list
#'
#' An object that stores a list of [cts_output] objects. Used when multiple
#' outputs are desired for a single module.
#'
#' @param ... Any number of [cts_output] objects.
#' @return An object that inherits from the class `cts_output_list`.
#'
#' @aliases cts_output_list
#' @export
output_list <- function(...) {
    structure(list(...), class="cts_output_list")
}

#' @rdname cts_output
#' @export
csv_output <- function(x, ..., name) {
    structure(x, class=c("cts_csv_output", "cts_output", class(x)), name=name)
}

#' @rdname cts_output
#' @export
graphical_output <- function(x, ..., name, width, height) {
    structure(x, class=c("cts_graphical_output", "cts_output", class(x)), name=name, width=width, height=height)
}

# #' @export
# tabular_output <- function(x, ..., name) {
#     structure(x, class=c("cts_tabular_output", "cts_output", class(x)), name=name)
# }
# 
# #' @export
# yaml_output <- function(x, ..., name) {
#     structure(x, class=c("cts_yaml_output", "cts_output", class(x)), name=name)
# }
# 
# #' @export
# rds_output <- function(x, ..., name) {
#     structure(x, class=c("cts_rds_output", "cts_output", class(x)), name=name)
# }

# @internal
cts_save <- function(x, ...) UseMethod("cts_save")

# @internal
#' @export
cts_save.cts_output_list <- function(x, ..., .path) {
    if (is.null(names(x))) {
        names(x) <- paste0("output", seq_along(x))
    }
    for (i in seq_along(x)) {
        set_name(x[[i]]) <- get_name(x[[i]]) %||% names(x)[i]
    }
    lapply(x, cts_save, .path=.path)
}

# @internal
#' @export
cts_save.cts_graphical_output <- function(x, ..., .path) {
    ggplot2::ggsave(plot=x, filename=file.path(.path, paste0(get_name(x), ".png")),
        width=(x %@% "width"), height=(x %@% "height"), ...)
}

# @internal
#' @export
cts_save.cts_csv_output <- function(x, ..., .path) {
    data.table::fwrite(x=x, file=file.path(.path, paste0(get_name(x), ".csv")), ...)
}

# # @internal
# cts_save.cts_tabular_output <- function(x, ..., .path) {
# }
# 
# # @internal
# cts_save.cts_yaml_output <- function(x, ..., .path) {
# }
# 
# # @internal
# cts_save.cts_rds_output <- function(x, ..., .path) {
# }



#********************************************************************************
#
#     UTILS
#
#********************************************************************************



# Functions copied from rlang
`%||%`  <- rlang::`%||%`
`%@%`   <- rlang::`%@%`
`%@%<-` <- rlang::`%@%<-`

merge_lists <- function(a, b) {

    if (is.null(b)) {
        return(a) 
    }

    # Extend the sequence x with m+1, m+2, ... so that it's length equals n
    ext_seq <- function(x, n, m) {
        lx <- length(x)
        if (n <= lx) {
            x[seq_len(n)]
        } else {
            c(x, 1 + seq.int(m, length.out=(n-lx)))
        }
    }

    nma <- names(a) %||% rep("", length(a))
    nmb <- names(b) %||% rep("", length(b))

    # Merge by position
    ia <- which(nma == "")
    ib <- which(nmb == "")
    ib <- ext_seq(ib, length(ia), length(b))
    b[ib] <- a[ia]

    # Merge by name
    nm <- setdiff(nma, "")
    b[nm] <- a[nm]

    b
}

merge_modules <- function(a, b) { 
    if (is.null(a) && is.null(b)) {
        return(NULL)
    }
    m <- merge_lists(a, b)
    for (x in c("module", "seed")) {
        attr(m, x) <- attr(a, x, exact=TRUE) %||% attr(b, x, exact=TRUE)
    }
    m
}

make_seed <- function() {
    round(1e8 * runif(1))
}



#' Counters
#'
#' Counters implemented as closures. Each time `next_counter()` is called,
#' it increments the counter.
#'
#' @param counter A character giving the name of the counter.
#' @param value Any value or object.
#' @param .all A logical. Should all counters be included?
#' @return `next_counter()` returns an integer.
#' @examples
#' next_counter("mycounter")
#' next_counter("mycounter")
#' next_counter("mycounter")
#' list_counters()
#' set_counter("mycounter", 10)
#' reset_counter("mycounter")
#' @name counter-functions
NULL

# @internal
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

#' @rdname counter-functions
#' @export
next_counter  <- closure_counter$next_counter

#' @rdname counter-functions
#' @export
set_counter   <- closure_counter$set_counter

#' @rdname counter-functions
#' @export
reset_counter <- closure_counter$reset_counter

#' @rdname counter-functions
#' @export
list_counters <- closure_counter$list_counters


#' Quantile interval
#'
#' A level \eqn{\alpha} quantile interval of `x` is the interval whose endpoints are the
#' upper and lower \eqn{(1-\alpha)/2} quantiles of `x`. In other words, it covers a
#' \eqn{(1-\alpha)} fraction of the data, with an equal number of points outside the
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

#' @import data.table
#' @importFrom stats quantile runif setNames update
#' @importFrom utils write.table
NULL

