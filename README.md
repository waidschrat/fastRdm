fastRdm - fast diffusion model analysis using R
====================================================

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

FastRdm estimates values for this parameters from response-time
distributions for both alternative decisions. FastRdm is based on
a fast numerical approach for solving the Partial Differential Equation 
defining the above sketched diffusion model (cf. Voss & Voss, 2008).
A detailed description of the core C functions that fastRdm relies
on is given by Voss and Voss (2007).

FastRdm comes with NO WARRANTY, to the extent permitted by law.  You
may redistribute copies of fastRdm under the terms of the GNU General
Public License.


References
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


Internal C Functions (fast-dm 30.2)
-----------------------------------

The core of the package is split into several source files. The common header
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


Users' manual
=============

Fast-dm is a program to efficiently estimate parameters in Ratcliff's
diffusion model.  Valid parameters are

    a = the boundary separation
    zr = the mean of the starting point of the diffusion relative to
         threshold separation
    szr = the width of the support of the distribution of zr
    v = the mean of the drift
    sv = the standard deviation of the drift
    t0 = the mean of the non-decisional component of the reaction time
    st0 = the width of the support of the distribution of t0
    d = differences in the non-decisional component between upper and
        lower threshold

The relative starting point is assumend to be uniformly distributed on
the interval [zr-0.5*szr,zr+0.5*szr], the drift is assumed to be
normally distributed with mean v and variance sv^2, and t0 is assumed
to be uniformly distributed on the interval [t0-0.5*st0,t0+0.5*st0].

For each run of the fast-dm program you need two kinds of input files:
one control file which describes the experiment and one or more data
files which contain the measured reaction times and responses.  In
addition there are two types of output files.  All files are plain
text files, the exact format is described below.


Experiment control files
------------------------

An experiment control file describes one experiment.  By convention
these files use the file name extension ".ctl".  By default fast-dm
tries to read an experiment description from the file
"experiment.ctl".  On systems which support command line arguments,
you can choose a different experiment control file by calling fast-dm
with the name of a control file as the first argument.

Control files are evaluated line by line.  Empty lines and lines
starting with "#" are ignored.  For all other lines the first word is
interpreted as a command and the subsequent words are used as
arguments to this command.

Example 1: A simple control file could look as follows

    format RESPONSE TIME
    load "*.dat"
    save "*.out"

The first line specifies that the data files will have just two
entries per line, namely first a 0/1 value giving the response and
then a reaction time (in seconds).  The second line instructs fast-dm
to analyse the data files "1.dat", "2.dat", etc.  The results are
written to the screen and to files named "1.out", "2.out", etc.

Example 2: A more complicated control file is the following one:

    precision 2
    set sv 0
    depends v stimulus
    format stimulus RESPONSE TIME
    load "*.dat"
    save "*.out"
    log "protocol"

Here, the first two lines instruct fast-dm to use reduced precision
for the computation and to not use variability for the v parameter (to
speed up the analysis).  The "format" statement specifies that each
line of the data files will contain three values, one entry describing
the stimulus used, then the participant's response and finally the
response time.  The "depends" statement specifies that different "v"
parameters will be used for different stimuli, all other parameters
are shared between all stimuli.  Datasets are read from the files
"1.dat", "2.dat", etc. and and the results are written to the files
"1.out", "2.out", etc. (in addition to being printed to the screen).
Also, a summary of all results, one line per input file, is written to
the file "protocol".


The following commands are valid in control files.  When present, they
should be given in the order of the following list.

  method M

    Use a given estimation method.  'M' can either be "ks" (for the
    Kolmogorov-Smirnov method; this is the default), "ml" (for the
    maximum likelihood method), or "cs" (for the chi-square method).

  precision VALUE

    Set the precision used for the computations in the parameter
    estimation procedure.  'VALUE' is a real number, the minimal
    allowed value is 1.  Larger values give higher accuracy but
    greatly increase computation time, reasonable values are in the
    range from 2 to 4.

    This command is optional, the default precision is 3.

  set PARAM VALUE

    Fixes the parameter 'PARAM' to 'VALUE' instead of estimating it.
    'PARAM' must be one of "a", "zr", "v", "t0", "d", "szr", "sv",
    st0".  Fixing parameters, especially the variability parameters,
    can speed up the estimation procedure significantly.

    This command is optional, default is to estimate all parameters.

  depends PARAM ...

    Specifies that parameter 'PARAM' depends on one or more
    experimental conditions.  'PARAM' must be one of "a", "zr", "v",
    "t0", "d", "szr", "sv", st0".  All the remaining arguments must be
    experimental conditions as defined on the "format" statement.

    Parameters which are fixed to a value using "set" cannot depend on
    experimental conditions.  All parameters which are not mentioned
    in a "depends" statement are shared between experimental
    conditions.

    This command is optional.  Default is to share all parameters
    between all experimental conditions.

  format ...

    Describes the format of the data files.  The arguments can be
    arbitrary strings.  Each argument corresponds to one column in the
    data files.  The arguments "RESPONSE" and "TIME" must occur
    exactly once each and define the position of the response column
    and reaction time column, respectively, in the datafile.  An
    argument of "*" specifies that the corresponding column is to be
    ignored.  All other arguments are interpreted as experimental
    conditions for use in "depends" statements.

    This command is required, at least RESPONSE and TIME must be
    specified.

  load "FNAME"

    Gives the name of the data files.  This command is required.  The
    given name may contain a single "*" to consecutively read several
    data files.

  save "FNAME"

    Gives the name of the per-dataset log files.  The results of
    parameter estimation, in addition to being printed to the screen,
    are written to these files.  If the file name in the 'load'
    contains a star, the file name given in the log command must also
    contain a star.

    This command is optional, default is not to write a per-dataset
    log file.  If this command is not present, "log" must be used
    instead.

  log "FNAME"

    Gives the name of the summary log file.  This command is optional,
    default is not to write a summary file.  If this command is not
    present, "save" must be used instead.



Data files
----------

Data files contain the participants' responses and measured reaction
times.  Optionally they can also contain information about
experimental conditions.  Typically each data file describes one
session of one participant in an experiment.  By convention data files
use the file name extension ".dat".

Data files are evaluated line by line.  Empty lines and lines starting
with "#" are ignored.  All other lines must contain the entries given
by the "format" statement in the control file.


Example 3.  The beginning of data file "17.dat" for "Example 2" above
could look as follows.

    # participant 17, female
    red 0 1.723
    black 1 0.877
    black 1 0.933
    red 1 1.526


The values in the "RESPONSE" column must be either 0 or 1,
corresponding to the two possible outcomes (a positive drift v or a
large z both favour the outcome 1).  The entries in the "TIME" column
must be floating point numbers, giving the reaction time in seconds.
All other entried correspond to experimental conditions and can be
arbitrary strings.



Per-Dataset Save Files
----------------------

per-dataset save files are written if the "save" statement is used in
the control file.  Each per-dataset log file corresponds to a data
file.  By convention, save files use the file name extension ".out".


Example 4.  The save file "17.out" for Example 3 above could look as
follows.

    a = 2.134450
    zr = 0.877730
    v_red = 2.52287
    v_black = -0.033281
    t0 = 0.322186
    d = 0.030334
    szr = 0.100000
    sv = 0.100000
    st0 = 0.100000
    p = 0.913227
    time = 42.850000

It gives the estimated parameters, one per line.  Since we allowed "v"
to depend on the stimulus, several values for "v" are given.  In
addition to the parameter values there are two informational entries:
'p' gives the resulting p-value (values close to one indicate a good
fit) and 'time' gives the CPU time used while processing this data set
in seconds.



Summary log file
----------------

A summary log file is created when the "log" command is used in the
control file.  This log file contains one line of output for each
dataset analysed plus a header line indicating which value is stored
in which column.  Summary log files can be useful to import fast-dm
results into other programs for statistical analysis.
