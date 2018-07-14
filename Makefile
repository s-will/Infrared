MAIN=infrared
#hello

HEADERS=assignment.hpp cluster.hpp cluster_tree.hpp		\
	constraint_network.hpp functions.hpp infrared.hpp	\
	rnadesign.hpp

CXXFLAGS=-Wall -pedantic

all: $(MAIN).so

check: test_$(MAIN)
	./test_$(MAIN)

run: all
	python3 ./test_$(MAIN).py

%.so: %.cpp $(HEADERS)
	g++ $(CXXFLAGS) --shared -fPIC -I /usr/include/python3.6 $< -lboost_python3 -o $@


test_$(MAIN): test_$(MAIN).cpp $(HEADERS)
	g++ $(CXXFLAGS) test_$(MAIN).cpp -o $@
