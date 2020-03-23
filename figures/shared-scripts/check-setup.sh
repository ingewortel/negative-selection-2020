
# Depending on system, check if basic tools are there.
system=$(uname -s)
if [ $system == "Darwin" ] ; then
	# This will open a prompt to install xcode CLT if it is not there, or otherwise
	# this will do nothing.
	xcode-select --install &> /dev/null
	echo "** 	Xcode CLT :	OK "
elif [ $system == "Linux" ] ; then
	checkEssentials=$( locate build-essential | wc -l )
	if [ $checkEssentials -eq 0 ] ; then
		echo "ERROR - Your computer does not appear to have build-essentials installed! Please install before continuing."
		echo ""
		echo "	Try:"
		echo "		sudo apt-get install build-essential"
		echo "	or try the equivalent for your package manager."
		exit 1
	else
		echo "** 	build-essentials :	OK "
	fi
else 
	echo "WARNING: Unknown system. Please ensure that you have basic command line tools (C compiler, make, ...) before continuing."
fi

# check latex
checkLatex=$(command -v pdflatex | wc -l )
checkLatexmk=$(command -v latexmk | wc -l )

if [[ $checkLatex == 0 || $checkLatexmk == 0 ]] ; then \
	echo "ERROR - Your computer does not appear to have latex (or all required packages) installed!"
	echo "	 Please install before continuing, and ensure you have the commands 'pdflatex', 'latexmk', and latex packages such as tikz installed."
	echo ""
	echo "	On linux, try:"
	echo "		sudo apt-get install texlive"
	echo "		sudo apt-get install texlive-latex-extra"
	echo "	On Mac OS X, try:"
	echo "		brew cask install mactex"
	exit 1
else
	echo "** 	Latex-etc :	OK "
fi

# check bc on linux (mac should have it)
checkBc=$(command -v bc | wc -l )

if [ $checkBc != 1 ] ; then \
	echo "ERROR - Your computer does not appear to have bc installed! Please install before continuing."
	echo ""
	echo "	On linux, try:"
	echo "		sudo apt-get install bc"
	echo "	On Mac OS X, you should have bc by default. Please Google to find out what's wrong."
	exit 1
else
	echo "** 	bc :		OK "
fi




# check openfst
checkOpenfst=$(command -v fstdifference | wc -l)

if [ $checkOpenfst != 1 ] ; then \
	echo "ERROR - Your computer does not appear to have OpenFST installed! Please install before continuing."
	echo ""
	echo "	On linux, try:"
	echo "		sudo apt-get install libfst-dev"
	echo "		sudo apt-get install libfst-tools"
	echo "	On Mac OS X, try:"
	echo "		brew install openfst"
	echo "	or visit http://www.openfst.org"
	exit 1
else
	echo "** 	OpenFST :	OK "
fi

# check graphviz
checkGraphviz=$(command -v neato | wc -l)

if [ $checkGraphviz != 1 ] ; then \
	echo "ERROR - Your computer does not appear to have Graphviz installed! Please install before continuing."
	echo ""
	echo "	On linux, try:"
	echo "		sudo apt-get install graphviz"
	echo "	On Mac OS X, try:"
	echo "		brew install graphviz"
	echo "	or visit https://www.graphviz.org/download/"
	exit 1
else
	echo "** 	Graphviz :	OK "
fi

# check rsvg-convert
checkRsvg=$(command -v rsvg-convert | wc -l)

if [ $checkRsvg != 1 ] ; then \
	echo "ERROR - Your computer does not appear to have rsvg-convert installed! Please install before continuing."
	echo ""
	echo "	On linux, try:"
	echo "		sudo apt-get install librsvg2-bin"
	echo "	On Mac OS X, try:"
	echo "		brew install librsvg"
	echo "	or Google how to install librsvg on your system."
	exit 1
else
	echo "** 	librsvg :	OK "
fi



# Check if R is installed
Rinstall=$(command -v R | wc -l)

if [ $Rinstall != 1 ] ; then \
	echo "ERROR - Your computer does not appear to have R installed! Please install R before continuing."
	echo ""
	echo "	On linux, try:"
	echo "		sudo apt-get install R"
	echo "	On Mac OS X, try:"
	echo "		brew install R"
	echo "	or visit https://cloud.r-project.org/"
	exit 1
else
	echo "** 	R :		OK "
fi


# Check if python 3 is installed
pythonCheck=$(command -v python3 | wc -l )

if [ $pythonCheck != 1 ] ; then \
	echo "ERROR - Your computer does not appear to have python3 installed! Please install python3 before continuing."
	echo ""
	echo "	On linux, try:"
	echo "		sudo apt-get install python3.6"
	echo "	On Mac OS X, try:"
	echo "		brew install python3"
	echo "	or visit https://www.python.org/downloads/"
	exit 1
else
	echo "** 	python3 :	OK "
fi



# install python packages
checkPythonPackages=$(python3 -m pip list)
checkNumpy=$(echo $checkPythonPackages | grep numpy | wc -l)

if [ $checkNumpy != 1 ] ; then \
	echo "You don't have the python3 Numpy packages. Do you wish to install it?"
	select yn in "Yes" "No"; do
		case $yn in
			Yes ) python3 -m pip install numpy; break;;
			No ) echo "ERROR - Please install numpy before continuing." ; exit;;
		esac
	done
else 
	echo "** 	 - Numpy :	OK "
fi

checkNetworkx=$(echo $checkPythonPackages | grep networkx | wc -l)

if [ $checkNetworkx != 1 ] ; then \
	echo "You don't have the python3 networkx packages. Do you wish to install it?"
	select yn in "Yes" "No"; do
		case $yn in
			Yes ) python3 -m pip install networkx; break;;
			No ) echo "ERROR - Please install networkx before continuing." ;exit;;
		esac
	done
else 
	echo "** 	 - Networkx :	OK "
fi



echo "Setup OK!"