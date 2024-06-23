#' conRes dataset
#'
#' This data sets contains the result that run from `BayesT` function using epil1 BCC object.
#' The epil1 object was obtained using `BCC.multi` function
#'
#' @usage data(conRes)
#' @format This is a dataframe with two columns and twenty observations
#' @examples
#' data(conRes)
#' conRes
"conRes"

#' PBCseqfit model
#'
#' This model contains the result that run from `BCC.multi` function using
#' PBC910 dataset in `mixAK` package
#'
#' @usage data(PBCseqfit)
#' @format This is a BCC model with thirty elements
#' @examples
#' data(PBCseqfit)
#' PBCseqfit
"PBCseqfit"

#' epil dataset
#'
#' This is epileptic.qol data set from `joinrRML`
#'
#' @usage data(epil)
#' @format This is a dataframe with 4 varaibles and 1852 observations
#' @examples
#' data(epil)
#' epil
"epil"

#' epil1 model
#'
#' This model contains the result that run from `BCC.multi` function using
#' epileptic.qol dataset in `joinrRML` package.
#' This model has formula of `formula =list(y ~ time + (1|id))`
#'
#' @usage data(epil1)
#' @format This is a BCC model with thirty elements
#' @examples
#' data(epil1)
#' epil1
"epil1"


#' epil2 model
#'
#' This model contains the result that run from `BCC.multi` function using
#' epileptic.qol dataset in `joinrRML` package.
#' This model has formula of `formula =list(y ~ time + (1 + time|id))`
#'
#' @usage data(epil2)
#' @format This is a BCC model with thirty elements
#' @examples
#' data(epil2)
#' epil2
"epil2"


#' epil3 model
#'
#' This model contains the result that run from `BCC.multi` function using
#' epileptic.qol dataset in `joinrRML` package.
#' This model has formula of `formula =list(y ~ time + time2 + (1 + time|id))`
#'
#' @usage data(epil3)
#' @format This is a BCC model with thirty elements
#' @examples
#' data(epil3)
#' epil3
"epil3"

#' example model
#'
#' This is an example model which contains the result that run from `BCC.multi`
#' function using epileptic.qol dataset in `joinrRML` package.
#' Only used in documented example and tests. Since small number of iterations
#' were used, this model can may not represent the true performance
#' for this method.
#'
#' @usage data(example)
#' @format This is a BCC model with thirty elements
#' @examples
#' data(example)
#' example
"example"

#' example1 model
#'
#' This is an example model which contains the result that run from `BCC.multi`
#' function using epileptic.qol dataset in `joinrRML` package.
#' Only used the tests. Since small number of iterations
#' were used, this model can may not represent the true performance
#' for this method.
#'
#' @usage data(example1)
#' @format This is a BCC model with thirty elements
#' @examples
#' data(example1)
#' example1
"example1"
