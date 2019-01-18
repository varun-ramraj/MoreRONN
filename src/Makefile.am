bin_PROGRAMS = moreRONN

#MoreRONN application

moreRONN_SOURCES = MoreRONN.cpp mclBBF.cpp
moreRONN_CPPFLAGS = -I./include -DDATA_PATH=\"$(prefix)\" -fopenmp
moreRONN_LDADD = -lrt -lm -lgomp

install-exec-hook:
	ln -fs $(prefix)/bin/moreRONN /usr/local/bin/moreRONN_49
	cp ../data/* $(prefix)

uninstall-hook:
	rm -rf $(prefix)
	rm -f /usr/local/bin/moreRONN_49
