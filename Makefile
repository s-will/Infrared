MAIN=sample
#hello

HEADERS=$(MAIN).hpp

all: $(MAIN).so

check: test_$(MAIN)
	./test_$(MAIN)

run: all
	./$(MAIN).py

%.so: %.cpp $(HEADERS)
	g++  --shared -fPIC -I /usr/include/python3.6 $< -lboost_python3 -o $@


test_$(MAIN): test_$(MAIN).cpp $(MAIN).hpp
	g++ test_$(MAIN).cpp -o $@
