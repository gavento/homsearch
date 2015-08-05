PYTHON=python
.PHONY: all clean test

all: homsearch_interface.so

clean:
	rm -rf homsearch_test build/ homsearch_interface.cpp homsearch_interface.so
	rm -rf homsearch.pyc homsearch_pytest.pyc

test: homsearch_test homsearch_interface.so
	./homsearch_test
	$(PYTHON) homsearch_pytest.py

homsearch_test: homsearch_test.cpp homsearch_lib.cpp homsearch_lib.h
	gcc -std=c++11 homsearch_test.cpp homsearch_lib.cpp -o homsearch_test -lstdc++ -g -Wall

homsearch_interface.so:	homsearch_lib.cpp homsearch_lib.h homsearch_interface.pyx
	$(PYTHON) setup.py build_ext --inplace
	rm -f homsearch_interface.cpp


