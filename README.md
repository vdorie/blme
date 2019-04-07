blme
====

Bayesian Linear Mixed Effect Models

A package for R. Built off of ['lme4'](https://cran.r-project.org/package=lme4)

Pre-built binaries of the package are available on [CRAN](https://cran.r-project.org/package=blme). These can be installed from within R using the typical `install.packages()` mechanism.

Steps to install from source:

1. Install development tools for your operating system (not necessary if you already have `lme4` installed):
    1. Linux/Unix should already have this installed; if not, use your package manager to install a C/C++ compiler.
    2. OS X: [clang and gfortran](https://cran.r-project.org/bin/macosx/tools/)
    3. Windows: [Rtools](http://cran.r-project.org/bin/windows/Rtools/)

2. Install the devtools package from within R:

```R
install.packages("devtools")
```

3. Run:

```R
install_github("vdorie/blme")
```
