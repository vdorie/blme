blme
====

Bayesian Linear Mixed Effect Models

A package for R. Built off of 'lme4' (http://cran.r-project.org/web/packages/lme4/index.html)

Pre-built binaries of the package are available on http://cran.r-project.org/web/packages/blme/index.html. These can be installed from within R using the typical `install.packages()` mechanism.

Steps to install from source:

  1. Install development tools for your operating system (not necessary if you already have lme4 installed):

    1. Linux/Unix should already have this installed
    2. OS X:
        1. Xcode (https://developer.apple.com/xcode/downloads/)
        2. gfortran (https://gcc.gnu.org/wiki/GFortranBinaries#MacOS)
    3. Windows: Rtools (http://cran.r-project.org/bin/windows/Rtools/)

  2. Install the devtools package from within R:

    `install.packages("devtools")`

  3. Run:

    `install_github("vdorie/blme")`
