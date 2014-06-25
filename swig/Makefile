# SWIG
INTERFACES = pysglite.i
WRAPPERS = $(INTERFACES:.i=_wrap.cxx)
PROXIES = $(INTERFACES:.i=.py )

all:
	swig -c++ -python  -I../ pysglite.i
	python setup.py install

# Clean target
.PHONY : clean
clean:
	$(RM) *.so
	$(RM) *.o
	$(RM) *.pyc
	$(RM) $(WRAPPERS)
	$(RM) $(PROXIES)