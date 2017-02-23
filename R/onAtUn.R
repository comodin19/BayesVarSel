##########################################################################
## start-up and clean-up functions
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 3, April 2013.
##
## Copyright (C) 2013-present Anabel Forte and Gonzalo Garcia-Donato
##    
##########################################################################

.onAttach <- function(...) {
	
	date <- date()
	x <- regexpr("[0-9]{4}", date)
	this.year <- substr(date, x[1], x[1] + attr(x, "match.length") - 1)
	
	# echo output to screen
	packageStartupMessage("#########")
	packageStartupMessage("## NOTE: Since v1.7.0 former function 'BayesFactor' has been renamed as 'Btest'")
	packageStartupMessage("## Copyright (C) 2013-", this.year,
			" Anabel Forte and Gonzalo Garcia-Donato", sep="")
	packageStartupMessage("#########")
}

.onUnload <- function(libpath) {
	library.dynam.unload("BayesVarSel", libpath)
}