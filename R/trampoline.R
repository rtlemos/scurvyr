#' Auxiliary function for tail recursions, part 1
#'
trampoline = function(f, ...) {
  function(...) {
    ret = f(...)
    while (inherits(ret, "recursion")) {
      ret = eval(as.call(c(f, unclass(ret))))
    }
    ret
  }
}

#' Auxiliary function for tail recursions, part 2
#'
recur = function(...) {
  structure(list(...), class = "recursion")
}
