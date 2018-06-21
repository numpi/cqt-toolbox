  ***************************************************************************
  * All the software  contained in this library  is protected by copyright. *
  * Permission  to use, copy, modify, and  distribute this software for any *
  * purpose without fee is hereby granted, provided that this entire notice *
  * is included  in all copies  of any software which is or includes a copy *
  * or modification  of this software  and in all copies  of the supporting *
  * documentation for such software.                                        *
  ***************************************************************************
  * THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED *
  * WARRANTY. IN NO EVENT, NEITHER  THE AUTHORS, NOR THE PUBLISHER, NOR ANY *
  * MEMBER  OF THE EDITORIAL BOARD OF  THE JOURNAL  "NUMERICAL ALGORITHMS", *
  * NOR ITS EDITOR-IN-CHIEF, BE  LIABLE FOR ANY ERROR  IN THE SOFTWARE, ANY *
  * MISUSE  OF IT  OR ANY DAMAGE ARISING OUT OF ITS USE. THE ENTIRE RISK OF *
  * USING THE SOFTWARE LIES WITH THE PARTY DOING SO.                        *
  ***************************************************************************
  * ANY USE  OF THE SOFTWARE  CONSTITUTES  ACCEPTANCE  OF THE TERMS  OF THE *
  * ABOVE STATEMENT.                                                        *
  ***************************************************************************

   AUTHORS:

       Dario A. Bini
       University of Pisa, Italy
       E-mail: bini@dm.unipi.it

       Stefano Massei
       EPFL, Switzerland
       E-mail: stefano.massei@epfl.ch

       Leonardo Robol
       ISTI-CNR, Italy
       E-mail: leonardo.robol@isti.cnr.it

   REFERENCE:

    -  Quasi-Toeplitz matrix arithmetic: a MATLAB toolbox
       NUMERICAL ALGORITHMS, NN (2018), pp. XXX-XXX

   SOFTWARE REVISION DATE:

       V1.0, June 2018

   SOFTWARE LANGUAGE:

       MATLAB R2011b and later


================================================================================
SOFTWARE
================================================================================

This software provides a MATLAB toolbox that implements the arithmetic of 
semi-infinite Toeplitz matrices with a low-rank correction. The toolbox 
provides operator overloading for the standard matrix operations, such as 
addition, multiplication, inversion, matrix powers, as well as a selection
of other useful matrix functions, such as the matrix square root and 
expoenntial. An implementatio of cyclic reduction to solve certain quadratic
matrix equations is included. 

================================================================================
PACKAGE
================================================================================

The package contains the following folders:

@cqt             - Implementation of the CQT class. 
docs             - Documentation and guide.
private          - A few internal implementations.
tests            - Unit tests.

On the top-level, the following files are available: 

cqtgallery.m      - Generates a few example matrices. 
cqtguide.m        - Compiles and open the guide.
cqtoption.m       - Set/Get various options for the toolbox.
cr.m              - Simple implementation of cyclic reduction.

================================================================================
HOW TO INSTALL
================================================================================

To install cqt-toolbox, it is sufficient to add it to your MATLAB path
by running addpath /path/to/cqt-toolbox

================================================================================
TEST
================================================================================

To run the unit tests, move the the tests folder and call the script RunAllTests.m

================================================================================
BUG FIXES / UPDATES
================================================================================

The most recent version of this package can be found on Github, at the
address https://github.com/numpi/cqt-toolbox. 
