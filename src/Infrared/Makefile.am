lib_LTLIBRARIES = libinfrared.la libhtdwrap.la

## install python modules
bin_SCRIPTS=treedecomp.py rna_support.py

########################################
## infrared module

libinfrared_la_SOURCES = libinfrared.cpp assignment.hpp		\
	boost_python_aux.hpp cluster.hpp cluster_tree.hpp	\
	constraint_network.hpp functions.hpp infrared.hpp	\
	rnadesign.hpp

########################################
## htd wrapper module

libhtdwrap_la_CXXFLAGS = $(AM_CXXFLAGS) $(HTD_CPPFLAGS)
libhtdwrap_la_LDFLAGS = $(LDFLAGS) $(HTD_LDFLAGS) $(HTD_LIBS)

libhtdwrap_la_SOURCES = libhtdwrap.cpp
