##############################################################################
#    Copyright (c) 2012-2016 Russell V. Lenth                                #
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

# Namespace hooks for emmeans

# .onLoad = function(libname, pkgname) {
#     # Function to check for existence of a variable
#     # This will be in the base package of R > 3.3.0
#     if (!exists("hasName", envir = getNamespace("utils"), inherits = FALSE)) {
#         assign("hasName", function(x, name)
#                 match(name, names(x), nomatch = 0L) > 0L,
#             envir = getNamespace(pkgname))
#     }
# }

# Just define the function for now. Maybe in a year we require R >= 3.4
#   (we now require that)
## hasName = function(x, name)
##    match(name, names(x), nomatch = 0L) > 0L


