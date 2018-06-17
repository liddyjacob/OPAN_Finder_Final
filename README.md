# OPAN Finder

This program finds all odd primitive abundant numbers with d distict prime factors.
This project is an implementation of: <http://ideaexchange.uakron.edu/honors_research_projects/728/>

We are in early stages of development, check back in a few weeks
and this will be more functional. If you would like to run this
program anyway, then keep reading.

Here are the things you will need:

* [NTL 11.0.0](http://www.shoup.net/ntl/)
* [GMP(for NTL)](https://gmplib.org/)
* [primecount](https://github.com/kimwalisch/primecount)

To run, use the following commands:

```bash
mkdir build
cd build
cmake ../
make
cmake ../
make
./run -help
```
