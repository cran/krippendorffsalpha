
######################################################################
#
# zzz.R
#
# Written by John Hughes <drjphughesjr@gmail.com>.
#
# Last Modified 10/11/22
# Licensed under the GNU General Public License version 2 (June, 1991)
#
######################################################################

#' @importFrom utils packageDescription
#' @importFrom utils flush.console
#' @importFrom stats quantile
#' @importFrom stats influence
#' @importFrom stats confint
#' @importFrom stats density
#' @importFrom stats qt
#' @importFrom stats var
#' @importFrom graphics plot
#' @importFrom graphics abline
#' @importFrom graphics hist
#' @importFrom graphics lines

.onAttach = function(libname, pkgname)
{
    temp = packageDescription("krippendorffsalpha")
    msg = paste("\n", temp$Package, ": ", temp$Title, "\n", "Version ", temp$Version,
                " created on ", temp$Date, ".\n", sep = "")
    msg = paste(msg, "copyright (c) 2020-2022, John Hughes\n",
                sep = "")
    msg = paste(msg, 'For citation information, type citation("krippendorffsalpha").\n', sep = "")
    msg = paste(msg, 'Type help(package = krippendorffsalpha) to get started.\n', sep = "")
    packageStartupMessage(msg)
}

