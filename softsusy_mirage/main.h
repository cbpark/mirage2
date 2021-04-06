
/** \file main.h
   - Project:     SOFTSUSY 
   - Author:      Ben Allanach 
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/
   - Description: main calling program header file: includes documentation
*/

/** \mainpage Detailed SOFTSUSY Documentation 

    \section install Installation or downloads
    For installation instructions or a download, please go to the 
    <a href="http://http://hepforge.cedar.ac.uk/softsusy/">
    SOFTSUSY Homepage</a>

    \section manual Official manual
    If you use SOFTSUSY to write a paper, please cite 
    <a href="http://xxx.soton.ac.uk/abs/hep-ph/0104145">
    B.C. Allanach, Comput. Phys. Commun. 143 (2002) 305-331, hep-ph/0104145, 
    </a> which is the SOFTSUSY manual.
    
    \section documentation Documentation
    These web-pages contain the documentation of the latest SOFTSUSY code.
    There are class diagrams and cross-referenced links a la doxygen to help 
    you navigate.

    \section updates Official Updates
    Updates will be posted on the    
    <a href="http://hepforge.cedar.ac.uk/softsusy/">
    SOFTSUSY Homepage</a>, and the 
    <a href="http://xxx.soton.ac.uk/abs/hep-ph/0104145">manual</a>
    will also be updated.
 */

#include <iostream>
#include <mycomplex.h>
#include <def.h>
#include <linalg.h>
#include <lowe.h>
#include <rge.h>
#include <softsusy.h>
#include <softpars.h>
#include <susy.h>
#include <utils.h>
#include <numerics.h>
#include <rpvsoft.h>
#include <rgeroutines.h>
#if defined(__hpux) || defined(_IBMR2) || defined(__alpha) || defined(__apollo) || defined(__ultrix) || defined(__sgi) || defined(__sun) 
#define extname
#endif 

#ifdef extname
#define FCN fcn_
#else
#define FCN fcn
#endif

#include "cfortran.h"
#include "minuit.h"
#define Ncont 20
