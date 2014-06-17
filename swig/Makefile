# SWIG
INTERFACES = sparsepad.i
WRAPPERS = $(INTERFACES:.i=_wrap.cxx)
PROXIES = $(INTERFACES:.i=.py )

# Default target: build the tests
.PHONY : all
all: $(WRAPPERS) ../sparsepad/sparsegrid.hpp sparsegrid.cxx
	./setup.py build_ext -i

# Rule: %.i -> %_wrap.cxx
%_wrap.cxx: %.i %.h numpy.i eigen.i
	swig -c++ -python $<
%_wrap.cxx: %.i %1.h %2.h numpy.i eigen.i
	swig -c++ -python $<

# Clean target
.PHONY : clean
clean:
	$(RM) *.so
	$(RM) $(WRAPPERS)
	$(RM) $(PROXIES)