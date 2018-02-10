#!/bin/sh

rm -rf autom4te.cache
libtoolize --force --copy && aclocal && autoheader && automake --add-missing --force-missing --copy && autoconf

#f77 -o testsym src/testsym.f
#f77 -c src/matfit.f
