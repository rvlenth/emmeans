##############################################################################
#    Copyright (c) 2012-2022 Russell V. Lenth                                #
#                                                                            #
#    This file is part of the emmeans package for R (*emmeans*)              #
#                                                                            #
#    *emmeans* is free software: you can redistribute it and/or modify       #
#    it under the terms of the GNU General Public License as published by    #
#    the Free Software Foundation, either version 2 of the License, or       #
#    (at your option) any later version.                                     #
#                                                                            #
#    *emmeans* is distributed in the hope that it will be useful,            #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of          #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           #
#    GNU General Public License for more details.                            #
#                                                                            #
#    You should have received a copy of the GNU General Public License       #
#    along with R and *emmeans*.  If not, see                                #
#    <https://www.r-project.org/Licenses/> and/or                            #
#    <http://www.gnu.org/licenses/>.                                         #
##############################################################################

# Extra examples

#' Run or list additional examples
#' 
#' This function exists so as to provide cleaner-looking examples in
#' help files when it must be run conditionally on another package.
#' Typically we want to run the code (\code{run = TRUE} is the default),
#' or otherwise just list it on the console (\code{list = TRUE}).
#'
#' @param name Character name of file to run. We look for a file with this name
#'   (with \code{".R"} appended) in the system files provided with \pkg{emmeans}.
#' @param run Logical choosing whether or not to run the example code
#' @param list Logical choosing whether or not to list the example code
#' @param ... Used only by the developer
#' 
#' @export
#' 
#' @examples
#' # List an example
#' emm_example("qdrg-biglm", list = TRUE)
#' 
#' # Run an example
#' if (require(biglm))
#'     emm_example("qdrg-biglm")
#' 
emm_example = function(name, run = !list, list = FALSE, ...) {
    test = (\(test = FALSE, ...) test)(...)
    if(test) {
        file = paste0("inst/extexamples/", name, ".R")
        filestg = paste0("inst/extexamples/", name, ".R")
    }
    else {
        file = system.file("extexamples", paste0(name, ".R"), package = "emmeans")
        if (file == "") stop("File '", name, ".R' not found")
        filestg = paste0("system.file(\"extexamples\", \"", paste0(name, ".R"),
        "\", package = \"emmeans\")")
    }
    if (list) {
        cat(readLines(file), sep = "\n")
        cat("\n")
    }
    if(run) {
        prompt.echo = try(get("prompt.echo", parent.frame(3)), silent = TRUE)
        continue.echo = try(get("continue.echo", parent.frame(3)), silent = TRUE)
        if (inherits(prompt.echo, "try-error")) {
            prompt.echo = getOption("prompt")
            continue.echo = getOption("continue")
        }
        message("\n--- Running code from '", filestg, "'")
        source(file, echo = TRUE, verbose = FALSE, 
               prompt.echo = prompt.echo, continue.echo = continue.echo)
        cat("\n")
    }
    invisible()
}
