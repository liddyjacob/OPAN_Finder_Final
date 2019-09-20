# OPAN Finder

This program finds all odd primitive abundant numbers with d distict prime factors.
This project is an implementation of: <http://ideaexchange.uakron.edu/honors_research_projects/728/>

If you would like to run this
program, you will need the following:

* [NTL 11.0.0](http://www.shoup.net/ntl/)
* [GMP(for NTL)](https://gmplib.org/)
* [primecount](https://github.com/kimwalisch/primecount)
* C++11
* cmake

To build and run the algorithm, open your terminal in 
this directory and use the following commands:

```bash
mkdir build
cd build
cmake ../
make
cmake ../
make
./run -help
```

TODO LIST:
* Build script in Python.
