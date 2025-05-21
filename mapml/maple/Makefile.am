
EXTRA_DIST = mapml.mpl mapml.mpl.bak maple-path.sh 

include_pmldir= $(includedir)/pml

# GV how to make mapml.mpl depend on mapml.mpl.bak ? Did not succeed for that 

install-data-local:
	sh maple-path.sh "$(prefix)"  "$(prefix)"
	$(INSTALL_DATA) mapml.mpl $(includedir)/pml/mapml.mpl
	


