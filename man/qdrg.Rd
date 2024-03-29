% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qdrg.R
\name{qdrg}
\alias{qdrg}
\title{Quick and dirty reference grid}
\usage{
qdrg(formula, data, coef, vcov, df, mcmc, object, subset, weights, contrasts,
  link, qr, ordinal, ...)
}
\arguments{
\item{formula}{Formula for the fixed effects}

\item{data}{Dataset containing the variables in the model}

\item{coef}{Fixed-effect regression coefficients (must conform to formula)}

\item{vcov}{Variance-covariance matrix of the fixed effects}

\item{df}{Error degrees of freedom}

\item{mcmc}{Posterior sample of fixed-effect coefficients}

\item{object}{Optional model object. \emph{This rarely works!}; 
but if provided, we try to set 
other arguments based on an expectation that `object` has a similar
structure to `lm` objects. See Details.}

\item{subset}{Subset of \code{data} used in fitting the model}

\item{weights}{Weights used in fitting the model}

\item{contrasts}{List of contrasts specified in fitting the model}

\item{link}{Link function (character or list) used, if a generalized linear model.
(Note: response transformations are auto-detected from \code{formula})}

\item{qr}{QR decomposition of the model matrix; used only if there are \code{NA}s
in \code{coef}.}

\item{ordinal}{\code{list} with elements \code{dim} and \code{mode}.
\code{ordinal$dim} (integer) is the number of levels in an ordinal response. If 
\code{ordinal} is provided, the intercept terms are modified appropriate to predicting 
an ordinal response, as described in \code{vignette("models")}, Section O,
using \code{ordinal$mode} as the \code{mode} argument (if not
provided, \code{"latent"} is assumed).
(All modes are supported except `scale`)
For this to work, we expect
the first \code{ordinal$dim - 1} elements of \code{coef} to be the
estimated threshold parameters, followed by the coefficients for the
linear predictor.}

\item{...}{Optional arguments passed to \code{\link{ref_grid}}}
}
\value{
An \code{emmGrid} object constructed from the arguments
}
\description{
This function may make it possible to compute a reference grid for a model 
object that is otherwise not supported.
}
\details{
Usually, you need to provide either \code{object}; or
\code{formula}, \code{coef}, \code{vcov}, \code{data}, and perhaps other
parameters. It is usually fairly straightforward to figure out how to get
these from the model \code{object}; see the documentation for the model class that
was fitted. Sometimes one or more of these quantities contains extra parameters,
and if so, you may need to subset them to make everything conformable. For a given \code{formula} and \code{data},
you can find out what is needed via \code{colnames(model.matrix(formula, data))}.
(However, for an ordinal model, we expect the first \code{ordinal.dim - 1} coefficients
to replace \code{(Intercept)}. And for a multivariate model, we expect \code{coef} 
to be a matrix with these row names, and \code{vcov} to have as many rows and columns as
the total number of elements of \code{coef}.)

If your model object follows fairly closely the conventions of an \code{\link[stats]{lm}}
or \code{\link[stats]{glm}}object, you may be able to get by providing the model as \code{object},
and perhaps some other parameters to override the defaults.
When \code{object} is specified, it is used as detailed below to try to obtain the 
other arguments. The user should ensure that the defaults
shown below do indeed work. 
The default values for the arguments are as follows:
\itemize{
  \item{\code{formula}: \code{formula(object)}}
  \item{\code{data}: \code{recover_data.lm(object)} is tried, and if an error is thrown,
    we also check \code{object$data}.}
  \item{\code{coef}: \code{coef(object)}}
  \item{\code{vcov}: \code{vcov(object)}}
  \item{\code{df}: Set to \code{Inf} if not available in \code{df.residual(object)}}
  \item{\code{mcmc}: \code{object$sample}}
  \item{\code{subset}: \code{NULL} (so that all observations in \code{data} are used)}
  \item{\code{contrasts}: \code{object$contrasts}}
}

The functions \code{\link{qdrg}} and \code{emmobj} are close cousins, in that
they both produce \code{emmGrid} objects. When starting with summary
statistics for an existing grid, \code{emmobj} is more useful, while
\code{qdrg} is more useful when starting from a fitted model.
}
\note{
For backwards compatibility, an argument \code{ordinal.dim} is invisibly 
supported as part of \code{...}, and if present, sets 
\code{ordinal = list(dim = ordinal.dim, mode = "latent")}
}
\section{Rank deficiencies}{

Different model-fitting packages take different approaches when the model
matrix is singular, but \code{qdrg} tries to reconcile them by comparing the
linear functions created by \code{formula} to \code{coefs} and \code{vcov}.
We may then use the \pkg{estimability} package to determine what quantities
are estimable. For reconciling to work properly, \code{coef} should be named
and \code{vcov} should have dimnames. To disable this name-matching
action, remove the names from \code{coef}, e.g., by calling \code{unname()}.
No reconciliation is attempted in multivariate-response cases. For more
details on estimability, see the documentation in the \pkg{estimability}
package.
}

\examples{
# In these examples, use emm_example(..., list = TRUE) # to see just the code

if (require(biglm, quietly = TRUE)) 
    emm_example("qdrg-biglm")
    
if(require(coda, quietly = TRUE) && require(lme4, quietly = TRUE)) 
    emm_example("qdrg-coda")
    
if(require(ordinal, quietly = TRUE)) 
    emm_example("qdrg-ordinal")

}
\seealso{
\code{\link{emmobj}} for an alternative way to construct an \code{emmGrid}.
}
