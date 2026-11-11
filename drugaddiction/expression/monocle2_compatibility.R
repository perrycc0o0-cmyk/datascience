# Compatibility with igraph 2.2.1 and scalar matrix checks in Monocle 2.
# These replacements are process-local and do not alter installed package files.
edit_call <- function(x, old, new) {
  if (is.call(x)) {
    if (is.symbol(x[[1]]) && identical(as.character(x[[1]]), old)) x[[1]] <- as.name(new)
    for (i in seq_along(x)) x[i] <- list(edit_call(x[[i]], old, new))
  }
  x
}

edit_argument <- function(x, old, new) {
  if (is.call(x)) {
    if (!is.null(names(x))) names(x)[names(x) == old] <- new
    for (i in seq_along(x)) x[i] <- list(edit_argument(x[[i]], old, new))
  }
  x
}

edit_projection_check <- function(x) {
  if (is.call(x)) {
    if (identical(x, quote(class(projection) != "matrix"))) return(quote(!is.matrix(projection)))
    for (i in seq_along(x)) x[i] <- list(edit_projection_check(x[[i]]))
  }
  x
}

apply_monocle2_compatibility <- function() {
  toy <- quote({ projection <- NULL; zz <- a[, 1]; z <- nei("a"); dfs(g, neimode = "all")
                 if (class(projection) != "matrix") projection <- as.matrix(projection) })
  expected <- quote({ projection <- NULL; zz <- a[, 1]; z <- .nei("a"); dfs(g, mode = "all")
                      if (!is.matrix(projection)) projection <- as.matrix(projection) })
  stopifnot(identical(edit_projection_check(edit_call(edit_argument(toy, "neimode", "mode"),
                                                       "nei", ".nei")), expected))
  ns <- asNamespace("monocle")
  update <- function(name, edit) {
    fn <- get(name, envir = ns, inherits = FALSE)
    body(fn) <- edit(body(fn))
    locked <- bindingIsLocked(name, ns)
    if (locked) unlockBinding(name, ns)
    on.exit(if (locked) lockBinding(name, ns), add = TRUE)
    assign(name, fn, envir = ns)
  }
  update("extract_ddrtree_ordering", function(x) edit_argument(x, "neimode", "mode"))
  update("project2MST", function(x) edit_projection_check(edit_call(x, "nei", ".nei")))
}
