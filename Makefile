MAIN=infrared
#hello

HEADERS=assignment.hpp cluster.hpp cluster_tree.hpp		\
	constraint_network.hpp functions.hpp infrared.hpp

all: $(MAIN).so

check: test_$(MAIN)
	./test_$(MAIN)

run: all
	./$(MAIN).py

%.so: %.cpp $(HEADERS)
	g++  --shared -fPIC -I /usr/include/python3.6 $< -lboost_python3 -o $@


test_$(MAIN): test_$(MAIN).cpp $(HEADERS)
	g++ test_$(MAIN).cpp -o $@
