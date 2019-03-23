Readme for DarkStars
p.scott@imperial.ac.uk

Short version of the standard licence:

  Copyright (c) 2009 Pat Scott
  Some residual parts of EZ Copyright (c) 2005, 2006,
  2007 Bill Paxton

  Permission is hereby granted, free of charge, to any
  person obtaining a copy of this software and associated
  documentation files (the "Software"), to deal in the
  Software without restriction, including without
  limitation the rights to use, copy, modify, merge,
  publish, distribute, sublicense, and/or sell copies of
  the Software, and to permit persons to whom the Software
  is furnished to do so, subject to the following
  conditions:

  The above copyright notice and this permission notice
  shall be included in all copies or substantial portions
  of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF
  ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED
  TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
  PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT
  SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
  CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
  OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR
  IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
  DEALINGS IN THE SOFTWARE.


Introduction
-------------

DarkStars is a stellar evolution package intended for use in
exploring the impact of dark matter upon stars.  It is based
on the EZ stellar evolution code, which is itself based on the
STARS code, and includes modified versions of some routines in
DarkSUSY.  The source code is freely available for download
from http://www.fysik.su.se/~pat/darkstars.

DarkStars has been described in a series of papers and conference
proceedings; if you publish work using or inspired by it, please
at least cite the two articles included in the distribution:
  Scott, Fairbairn & Edsjö (2009, MNRAS 394:82, arXiv:0809.1871)
  Scott, Edsjö & Fairbairn (2009, Proc. DARK2009, arXiv:0904.2395).
Please refer to those two articles for further information on the
DarkStars code and its input physics.  You may also consider citing
one or more of the other articles describing DarkStars:
  Fairbairn, Scott & Edsjö (2008, PRD 77:047301, arXiv:0710.3396)
  Scott, Edsjö & Fairbairn (2008, Proc. DARK2007, arXiv0711.0991)
  Scott, Fairbairn & Edsjö (2008, Proc. IDM08, PoS(idm2008)073,
                            arXiv:0810.5560).

If you publish something using DarkStars, please also cite the
EZ paper:
  Paxton (2004, PASP 115:699, arXiv:astro-ph/0405130),
as well as the original STARS papers:
  Eggleton (1971, MNRAS 151:351)
  Eggleton (1972, MNRAS 156:361)
  Pols, Tout, Eggleton & Han (1995, MNRAS 274:964, arXiv:astro-ph/9504025).


Installation
-------------


Unpack the tar file to a directory of your choice, cd to DarkStars/run
and type './mkdark'.  If you've got the Intel or GNU fortran compiler,
and are running on an x86 machine, then everything
*should* go smoothly.  If not, edit DarkStars/make/makefile
to suit your compiler and architecture.  Users should probably do this
anyway to ensure that the correct optimisations are being used for their
specific architecture.  It may also be necessary to edit the makefiles
in the DarkSUSY directory if compiling proves troublesome.

DarkStars has been tested with ifort 12 and gfortran 4.8.5.


Running (and plotting) the examples
-----------------------------------
DarkStars operates by piping simple text files into the program as input.
As output, it writes a series of log files, prints varying degrees of
information to STDOUT and sometimes dumps .sav and .sav.dark files.

To run the examples, try first

cd DarkStars/run
./DarkStars < ../DarkRuns/example_verydark/example_verydark.in > ../DarkRuns/example_verydark/example_verydark.out

This run should be reasonably quick, i.e. about 30 minutes.  For the other
two examples, replace 'verydark' in the command above with 'dark' or 'normal'.
Be aware that the 'dark' example is set to run at 10 times the temporal
resolution of the other two examples, so will probably run for over a day!

When you're done running examples, compare the output in the
DarkStars/DarkRuns/example_*/logs/*.log and DarkStars/DarkRuns/example_*/example_*.out
files to the example output in DarkStars/DarkRuns/example_*/diffme.

The DarkStars/DarkRuns/example_*/diffme/plots directory also contains example plots
for each example star.  If you have Ruby and Tioga properly installed you can
generate your own plots by editing DarkStars/star_profile/path.txt to point to
the run you're interested in, then going to DarkStars/star_history and running

ruby make_history_pdfs.rb

The plots will appear in DarkStars/star_history/history_out.  Other plots can
be generated using DarkStars/star_profile/make_profile_pdfs.rb and
DarkStars/PatPlots/make_Pat_pdfs.rb.  Have a poke around in the rb files in
these directories if you want to modify the plots (I seem to need to change
spacings, font sizes, legends etc all the time because I want a slightly different
version of a plot - so since some plots use common code, a few of them are not
currently spaced very nicely).

There are a range of options and inputs available to the user in a DarkStars input
file - the following sections give a quick explanation of what each one does.
Various of these are incompatible with each other, so the code checks for this and
throws a (hopefully informative) error message when a conflict occurs.  If any are
not specified, they are assumed to be zero - which is not checked for in the code,
so be careful not to leave any out!


Directories
-------------
 metals_dir: directory to use for metallicity information.
             Should be one of z0_proto, z0001, z0003, z001,
             z004, z01, z02, or z03

 main_dir: the main directory for I/O for this model.
           This is where the data, profile and movie
           directories must reside, and will be where
           any save/restore files are created/sourced.
           This is also where DarkStars will search for
           orbit.dat if the user chooses to specify a
           particular orbit.

 data_dir: the directory in which to save output log files.
           The format of the standard EZ logs can be found in
           DarkStars/demos/DATA_README.  The format of the
           additional log added by DarkStars (extra.log) can
           be found (as comments) and modified in lines 260-290
           of ez_do_one.f.

 prof_dir: the directory in which to save profile files
           for important individual times in the star's
           evolution.

 movie_dir: the directory in which to save profile files
            at regular increments for use in creating
            movies.


Switches
--------

 capmode
false: analytical integration over Gaussian
       halo profile
 true: numerical calculation

 interpmode
false: linear interpolation in mass, density
       fractions, potential and WIMP radial
       number density.  This mode has not been
       fully implemented yet.
 true: cubic spline interpolation in some things
       and tensional spline interpolation using
       TSPACK in others.

 potmode
false: tabulation of potentials, with
       interpolation as per interpmode
 true: explicit evaluation of potential
       integrals during integration over
       halo profile and annihilation rates,
       as well as evaluation of the WIMP
       gravitational partition function.
 NOTE: tabulation is faster, but might produce inaccurate
       results if the EZ adaptive mesh gets really
       wild, particularly if interpmode is set to true.

 chopatvesc
false: do not assume any cut-off in the halo velocity
       distribution.
 true: truncate the WIMP velocity distribution at the
       local galactic escape velocity.

 altveldist
false: assume a standard isothermal WIMP halo velocity
       distribution.
 true: use some other fancy user-defined WIMP velocity
       distribution.  Currently defaults to the 'N-body'
       distribution described in the 2009 MNRAS paper.

 WIMPdensmode
false: approximate expression used for radial
       WIMP distribution, based upon the
       assumption of constant stellar density
       in the region populated by WIMPs (i.e.
       the central density).
 true: WIMP radial distribution calculated
       based upon true density distribution.

 ReconvergeModels
false: run in explicit mode whilst solving capture
       and stellar structure equations.
 true: reconverge each stellar model after each
       timestep with the new WIMP population and
       distribution.

 DoTransport
false: WIMP conductive energy transport
       ignored
 true: WIMP energy transport enabled

 DoMovie
false: don't save files required for making an evolutionary
       movie
 true: save required data for creating a movie

 DoWIMPyThings
false: don't compute effects, capture nor annihilation of WIMPs
       at all
 true: go to town

 verbose
false: quiet output.  Errors, warnings and comments are still
       printed.
 true: verbose output, printing WIMP population, capture rates,
       etc at every timestep.

 Integrator
 0:    except for capture integrals, all integration
       is performed using Simpson's rule. This is the
       default.
 1:    use Romberg integration for everything except
       capture integrals.  This option has been
       disabled in the public release as it employs
       a proprietary integrator, but the user can
       easily re-enable it with their own version.
       Look in the file /src/DkStrs_fint_ext.f
       for more details.
 2:    use Runge-Kutta for everything except capture
       integrals.  Also disabled; see option 1.


Controls
--------

 hold_on:  the target maximum fractional change in
           the WIMP population per timestep. This is
           used to limit the timestep, but can be
           overridden at certain times, such as in
           the first few timesteps of any run.

 initial_timestep: the length of the first 3 timesteps
                   in the run, in years.  The code can
                   be *very* sensitive to the value of
                   this control, and it can often take
                   a lot of fine tuning of the initial
                   timestep to get a particular model
                   to converge correctly.

 timestep_rescale: the scaling factor to apply to the
                   internal STARS timestep chooser.
                   timestep_rescale = 5.5*CDD where CDD
                   is the control factor in the original
                   STARS implementation.

 min_timestep: the minimum timestep allowed when a star
               follows a user-defined orbit given by
               orbit.dat.  This prevents timesteps getting
               ridiculously small and ruining convergence
               when e.g. a star is passing through periapsis
               on a very tight elliptical orbit.

 summary_cnt: save log and movie data every summary_cnt
              timesteps.

 stop_at_Model: if positive, this is the model number to stop and
                save at for later use.  There are a number of
                negative integral values which constitute special
                stop codes:
                -1
                -10 Stop when the star either leaves the main
                    sequence or becomes a dark star (halts in
                    the HR diagram).  This is just a logical
                    AND of -6 and -x
                Any other negative value is interpreted as an
                indication that the code should not stop and
                save its state for later use at any particular
                model number.

 max_Age: an age to stop and save at for later use. Negative
          values are interpreted as an indication that the code
          should not stop and save its state for later use at
          any particular stellar age.  Should be given in
          years.

 save_filename: the filename to use in the case where the code
                will save its state for later use.

 restore_filename: the filename to restore a previously saved
                   model from.  Set this to '' for a new run.

 save_freq: The time in days to run between periodic saves. Set
            this negative to just save at the end (as specified
            by max_Age or stop_at_Model).


Physical Inputs
---------------

 mx: the WIMP mass in GeV.

 sigsi: the spin-independent WIMP-nucleus scattering cross-section
        in cm^2.

 sigsd: the spin-dependent WIMP-nucleus scattering cross-section
        in cm^2.

 sigann: the velocity-averaged WIMP self-annihilation cross-section,
         in the limit v->0, in cm^3/s.

 nuLossFactor: the final fraction of WIMP annihilation energy lost
               to neutrinos (dimensionless).

 boost_factor: a factor by which capture rates will be artificially
               multiplied in order to mimic some sort of ignored
               physical effects (dimensionless).

 vd_3d_star: the velocity dispersion of the dark matter halo in km/s.

 v_star: the stellar proper velocity through WIMP halo in km/s. Set
         this negative to obtain this quantity dynamically from
         orbit.dat.

 rhowimp: the local ambient WIMP density in GeV/cm^3. Set this
          negative to obtain this quantity dynamically from orbit.dat.

 galesc: the local galactic escape velocity in km/s. Set
         this negative to obtain this quantity dynamically from
         orbit.dat.

 galr: the local galactocentric distance of the star in parsecs.  Set
         this negative to obtain this quantity dynamically from
         orbit.dat.

 n_WIMPs: the initial population of WIMPs in the star with which to
          begin the simulation.

 mass: the initial stellar mass, in solar masses.  This can be
       anything from 0.3 to 100.

 K0: the value of the Knudsen number at which WIMP energy transport
     is maximised (dimensionless).

 Knudsen_suppression_tau: the relaxation scale for the Knudsen-
                          dependent suppression function (dimensionless)

 tau_therm: the WIMP thermalisation timescale in seconds. Set negative
            to assume that WIMPs automatically thermalise in the star
            as soon as they are captured.  This option is very
            'experimental'.  Setting it to a positive value forces DarkStars
            to use only the approximation for the WIMP distribution
            in the star specified by WIMPdensmode = false, ignore
            any corrections to the density structure due to a finite
            WIMP mean-free path, and ignore all WIMP conductive energy
            transport.

Abundances
-----------
These specify the solar elemental abundances and ratios used for the legacy
solar capture mode (accessible only via a logical switch in DkStrs_constants.f),
and for determining how the heavy-element abundances rescale with
metallicity.  Absolute abundances are given in logarithmic units log_10(x/H)+12,
where x and H are the number densities of element x and hydrogen respectively,
and the abundance of hydrogen is defined as 12.  Isotopic ratios are direct
ratios of number densities.  Most users will have no reason to modify these
values from those used in the example file.

 solar_C: the solar carbon abundance in logarithmic units.

 solar_O: the solar oxygen abundance in logarithmic units.

 heavy_abuns: the solar Na, Al, Si, S, Ar, Ca, Fe, Ni & Pb abundances in
              logarithmic units.

 iso_ratios: the solar ratios of He4/He3, C12/C13, O16/O18, Ni58/Ni60,
             Pb208/Pb207 & Pb208/Pb206.


Advanced options: user-specified orbits
-----------------------------------------
A user can provide an explicit orbit for their star to follow in a
file 'orbit.dat' placed in the main_dir directory.  Examples of these
files are given in DarkStars/orbitdata.  The orbit.dat file must begin
with an integer giving the total number of points at which the orbit is
sampled, followed by a real number specifying the minimum fraction of
the orbit over which the code is allowed to assume that the ratio of
the ambient WIMP density to the stellar velocity changes monotonically.
This is used as an upper limit for timesteps, in order to prevent the
code accidentally skipping over transient peaks or troughs in the
WIMP density or stellar velocity by choosing a timestep purely on the
basis of the current stellar structure.  These are then followed by
columns specifying the dynamical variables at each sampled point
in the orbit.  The columns are:
time(yr)
WIMP density(GeV/cm^2)
stellar velocity (km/s)
local Galactic escape velocity (km/s)
local Galactocentric distance (cm).

In order for the orbit.dat file to be used, at least rhowimp and
v_star must be set negative in the input file.  If galesc and galr
are set positive in the input file, their values in orbit.dat will
be ignored in favour of the constant values in the input file.
If one of these two quantities is specified as negative in the input
file, indicating that it should be determined dynamically from
orbit.dat, then the other must be as well - i.e. the only options are
all 4 positive
all 4 negative
rhowimp & v_star negative, galesc & galr positive.

If the final row and the first row in orbit.dat have exactly the same
values in all columns except time, then the orbit will loop.


Advanced options: Z = 0 protostellar evolution
-----------------------------------------------
In order to run a metal-free protostar evolutionary simulation, the user
needs to specify metals_dir = '../metals/z0_proto' and supply a
protostellar starting model.  In this case, the input 'mass' parameter
will be ignored.  Starting models can be input via
DarkStars/metals/z0_proto/protostar.start.  This file contains
199 lines corresponding to the meshpoints located at different heights in
the star, and currently contains a 20 solar mass ZAMS starting model as an
example (not a Z = 0 model though!).  Columns correspond to the values of
the following quantities at each meshpoint:
ln f (degeneracy parameter)
ln T(Kelvin)
the abundance (mass fraction) of 16O
ln m - mass interior to this meshpoint in units of 10^30kg
the abundance of 1H
the gradient of the mesh spacing function
ln r - the radius in units of 10^9m
L - luminosity in units of 10^26W
the abundance of 4He
the abundance of 12C
the abundance of 20Ne

As a related aside, the WIMP distribution inside the star is provided by
the routine WIMPdens in EZ/src/DkStrs_WIMPdens.f - this is what must be
altered if one wanted to say, allow an adiabatically-contracted density
profile instead of a thermalised one.


Advanced options: user-specified velocity distributions
-------------------------------------------------------
In order to specify their own dark matter velocity distribution, the user
need simply write their own function returning the value of their distribution
function for some velocity v, and modify the function alt_foverusq in
DarkStars/DarkSUSLE/src/cap/dscapfoveru.f to point to it.  This distribution
can then be called using altveldist = true, with the option also of
automatically truncating it at the local escape velocity with
chopatvesc = true.  Any required renormalisation of the distribution will be
taken care of by DarkStars in this case.

