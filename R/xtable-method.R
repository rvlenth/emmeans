##############################################################################
#    Copyright (c) 2012-2017 Russell V. Lenth                                #
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

### xtable method
# Modified from xtableLSMeans function provided by David Scott

#' Using \code{xtable} for EMMs
#' 
#' These methods provide support for the \pkg{xtable} package, enabling 
#' polished presentations of tabular output from \code{\link{emmeans}}
#' and other functions.
#' 
#' The methods actually use \code{\link[xtable]{xtableList}},
#' because of its ability to display messages such as those for P-value
#' adjustments. These methods return an object of class \code{"xtable_emm"} --
#' an extension of \code{"xtableList"}. Unlike other \code{xtable} methods, the
#' number of digits defaults to 4; and degrees of freedom and \emph{t} ratios
#' are always formatted independently of \code{digits}. The \code{print} method
#' uses \code{\link[xtable:xtableList]{print.xtableList}}, and any \code{\dots} arguments are
#' passed there.
#'
#' @param x Object of class \code{emmGrid}
#' @param caption Passed to \code{\link[xtable]{xtableList}}
#' @param label Passed to \code{xtableList}
#' @param align Passed to \code{xtableList}
#' @param digits Passed to \code{xtableList}
#' @param display Passed to \code{xtableList}
#' @param auto Passed to \code{xtableList}
#' @param ... Arguments passed to \code{\link{summary.emmGrid}}
#'
#' @return The \code{xtable} methods return an \code{xtable_emm}
#' object, for which its print method is \code{print.xtable_emm} .
#' 
#' @method xtable emmGrid
#' @importFrom xtable xtable
#' @importFrom xtable xtableList
#' @export
#' @examples
#' pigsint.lm <- lm(log(conc) ~ source * factor(percent), data = pigs)
#' pigsint.emm <- emmeans(pigsint.lm, ~ percent | source)
#' xtable::xtable(pigsint.emm, type = "response")
xtable.emmGrid = function(x, caption = NULL, label = NULL, align = NULL, digits = 4, 
    display = NULL, auto = FALSE, ...) 
{
    xtable.summary_emm(summary(x, ...), caption = caption, label = label, align = align, digits = digits, 
           display = display, auto = auto)
}

#' @rdname xtable.emmGrid
#' @method xtable summary_emm
#' @export
xtable.summary_emm = function (x, caption = NULL, label = NULL, 
          align = NULL, digits = 4, 
          display = NULL, auto = FALSE, ...) 
{
    if (!is.null(x$df)) x$df = round(x$df, 2)
    if (!is.null(x$t.ratio)) x$t.ratio = round(x$t.ratio, 3)
    if (!is.null(x$z.ratio)) x$z.ratio = round(x$z.ratio, 3)
    if (!is.null(x$p.value)) {
        fp = x$p.value = format(round(x$p.value,4), nsmall=4, sci=FALSE)
        x$p.value[fp=="0.0000"] = "<.0001"
    }
    if (!is.null(byv <- attr(x, "by.vars"))) {
        byc = which(names(x) %in% byv)
        xList = split(as.data.frame(x), f = x[, byc])
        labs = rep("", length(xList))
        for (i in 1:length(xList)) {
            levs = sapply(xList[[i]][1, byc], as.character)
            labs[i] = paste(paste(byv, levs, sep = " = "), collapse = ", ")
            xList[[i]] = as.data.frame(xList[[i]][, -byc, drop = FALSE])
        }
        attr(xList, "subheadings") = labs
    }
    else {
        xList = list(as.data.frame(x))
    }
    attr(xList, "message") = attr(x, "mesg")
    result = xtable::xtableList(xList, caption = caption, label = label, 
       align = align, digits = digits, display = display, 
       auto = auto, ...)
    digits = xtable::digits(result[[1]])
    
    # format df and t ratios
    digits = xtable::digits(result[[1]])
    i = which(names(x) == "df")
    if (length(i) > 0) {
        dfd = ifelse(all(zapsmall(x$df - round(x$df)) == 0), 0, 2)
        digits[i + 1 - length(byv)] = ifelse(is.na(dfd), 0, dfd)
    }
    i = which(names(x) %in% c("t.ratio", "z.ratio"))
    if (length(i) > 0) digits[i + 1 - length(byv)] = 3
    for (i in seq_along(result))
        xtable::digits(result[[i]]) = digits
    
    class(result) = c("xtable_emm", "xtableList")
    result
}

#' @rdname xtable.emmGrid
#' @param type Passed to \code{\link[xtable]{print.xtable}}
#' @param include.rownames Passed to \code{print.xtable}
#' @param sanitize.message.function Passed to \code{print.xtable}
#' @method print xtable_emm
#' @export
print.xtable_emm = function(x, type = getOption("xtable.type", "latex"),
                            include.rownames = FALSE, 
                            sanitize.message.function = footnotesize,
                            ...)
{
    footnotesize = switch(type,
        html = function(x) paste0("<font size = -1>", x, "</font>"),
        latex = function(x) paste0("{\\footnotesize ", x, "}"),
        function(x) x )
    invisible(xtable::print.xtableList(x, type = type, include.rownames = include.rownames, 
        sanitize.message.function = sanitize.message.function, ...))
}
