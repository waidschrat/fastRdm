fast-dm - a fast method for diffusion model analysis


INTRODUCTION
------------

Stochastic Diffusion Models are used in cognitive science to analyse
cognitive processes in fast binary decisions. It is assumed that the 
information supporting one and refuting the other decision can be
described by a Wiener Diffusion process with constant drift over time.
In the model, the decision process is terminated when the process 
exits from a given interval (see also Voss, Nagler, & Lerche, 2013).

The diffusion-model data analysis  was introduced in psychology by
Roger Ratcliff (1978). Ratcliff's diffusion model is specified by the
following parameters: The distance of thresholds (a), the starting
point (z), the drift rate (v), a non-decisional part of the response
time (t0), and three so-called inter-trial variability parameters
concerning the starting point (sz), the drift (sv), and the non-
decisional response-time constant (st0).

Fast-dm estimates values for this parameters from response-time
distributions for both alternative decisions. Fast-dm is based on
a fast numerical approach for solving the Partial Differential Equation 
defining the above sketched diffusion model (cf. Voss & Voss, 2008).
A detailed description of fast-dm is given by Voss and Voss (2007).

Instructions about how to use the programm can be found in the
file MANUAL.  Please mail any suggestions and bug reports to
Jochen Voss <voss@seehuhn.de> or to
Andreas Voss <Andreas.Voss@psychologie.uni-heidelberg.de> .


Fast-dm comes with NO WARRANTY, to the extent permitted by law.  You
may redistribute copies of fast-dm  under the terms of the GNU General
Public License.  For more information about these matters, read the
file COPYING of the source code distribution.


REFERENCES
----------
Ratcliff, R., 1978. A theory of memory retrieval. Psychological Review,
  85, 59-108.

Ratcliff, R., Rouder, J. N., 1998. Modelling response times for two-choice
  decisions. Psychological Science 9 (5), 347-356.

Voss, A., Rothermund, K., Voss, J., 2004. Interpreting the parameters of the
  diffusion model: An empirical validation. Memory & Cognition 32, 1206-1220.

Voss, Nagler, & Lerche (in press). Diffusion Models in Experimental Psychology:
  A Practical Introduction. Experimental Psychology.

Voss, A., Voss, J., 2008. A fast numerical algorithm for the estimation
  of diffusion-model parameters. Journal of Mathematiocal Psychology.

Voss, A., Voss, J., 2007. Fast-dm: A free program for efficient diffusion model
  analysis.  Behavioral Research Methods, 39, 767-775.

Wagenmakers, E.-J., van der Maas, H. L. J., & Grasman, R. P. P. P., 2007.
  An ez-diffusion model for response time and accuracy. 
  Psychonomic Bulletin & Review, 14, 3-22.


INSTALLATION
------------

On Unix-like systems:
Just type the following commands

    ./configure
    make

On Windows systems:
As usual, on a Windows system things are a little more complicated.  The
following steps may help you if you are using Microsoft's Visual Studio 7:

(1) Create a new project. Choose a "win32 Console project", and make sure 
    that the "empty project" check-box is selected.
(2) Copy all "fast-dm" files in the new project directory.
(3) Add the source code to the project.  You need the files EZ-diff.c,
    cdf.c, container.c, dataset.c, experiment.c, file.c, main.c,
    method-ks.c, pde.c, phi.c, simplex2.c, win32dir.c, win32dir.h,
    win32erf.c, xmalloc.c and fast-dm.h .
(4) In the project settings, select "compile as C-code" (and not C++).
(5) Compile the project.


INTERNALS
---------

The program is split into several source files.  The common header
file for all of these is "fast-dm.h".

    cdf.c         - compute the CDF for the diffusion model
    container.c   - container types (sets, dictionaries, arrays)
    dataset.c     - read and store data files
	density.c     - calculate the density of the predicted response-time distributions
    experiment.c  - read and store the control file
    EZ-diff.c     - compute approximate parameters using the EZ-model
    file.c        - auxiliary funcitons to read control and data files
    main.c        - main program for the fast-dm project
    method-cs.c   - implementation of the CS method
    method-ks.c   - implementation of the KS method
    method-ml.c   - implementation of the ML method
    pde.c         - numerically solve the Fokker-Planck equation
    phi.c         - the CDF and inverse CDF of the standard normal distribution
    precision.c   - tune parameters for the precision of CDF and density
    simplex2.c    - the downhill simplex method of Nelder and Mead
    win32dir.c    - emulate opendir/readdir/closedir on MS Windows systems
    win32.erf.c   - emulate Gaussian error function (erf) on MS Windows systems
    win32getopt.c - emulate getopt on MS Windows systems
    xmalloc.c     - memory management

    plot-cdf.c - plot cumulative distribution functions

    construct-samples.c - simulate sample data
