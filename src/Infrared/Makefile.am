lib_LTLIBRARIES = libinfrared.la

## install python modules
bin_SCRIPTS=treedecomp.py rna_support.py infrared.py

########################################
## infrared headers
infrared_includedir=$(includedir)/Infrared
infrared_include_HEADERS=assignment.hpp boost_python_aux.hpp	\
	cluster.hpp cluster_tree.hpp constraint_network.hpp	\
	functions.hpp infrared.hpp rnadesign.hpp

########################################
## infrared module

libinfrared_la_SOURCES = libinfrared.cpp $(infrared_include_HEADERS)

########################################
## htd wrapper module
if HAVE_HTD
lib_LTLIBRARIES+=libhtdwrap.la
libhtdwrap_la_CXXFLAGS = $(AM_CXXFLAGS) $(HTD_CPPFLAGS)
libhtdwrap_la_LDFLAGS = $(LDFLAGS) $(HTD_LDFLAGS) $(HTD_LIBS)

libhtdwrap_la_SOURCES = libhtdwrap.cpp
endif
