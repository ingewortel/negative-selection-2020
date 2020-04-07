
argv <- commandArgs( trailingOnly = TRUE )
Rpackfile <- argv[1]

# Read list of packages from file
package.list <- read.table( Rpackfile )$V1

cat("Checking computer for required R packages : \n")

not_installed <- character()
for( p in package.list ){

	# See if package is installed
	test <- suppressMessages( suppressWarnings( require( p, character.only = TRUE, quietly = TRUE ) ) )

	# print diagnostic info
	cat("....",p,"\t: ",test, "\n")

	# add to list of not installed if test = FALSE
	if( !test ){
		not_installed <- c( not_installed, p )
	}

}

if( length( not_installed ) == 0 ){

	message( "\n==> Setup OK. \n")

} else {

	message( paste( "\nThe following packages are missing on your computer: \n\t", not_installed,
		"\nDo you wish to install them now? y = yes, n = no." ) )

	a <- readLines("stdin", n=1)
	if( a == "y" ){
		for( p in not_installed ){
			message( "Installing package: ", p )
			suppressMessages( install.packages( p, quiet=TRUE ) )
		}		


		message( "\n\n\n==> Setup OK. \n")

	} else {
		stop( "Please install missing packages manually before continuing." )
	}

}
