J/A+AS/130/65       A standard stellar library. II. (Lejeune+ 1998)
================================================================================
A standard stellar library for evolutionary synthesis:
II. The M dwarf extension
       Lejeune T., Cuisinier F., Buser R.
      <Astron. Astrophys. Suppl. Ser. 130, 65 (1998)>
      =1998A&AS..130...65L      (SIMBAD/NED BibCode)
================================================================================
ADC_Keywords: Atlases ; Models, evolutionary ; Stars, late-type ;
              Photometry, UBVRIJKLMNH
Keywords: stars: fundamental parameters - stars: late type

Description:
    A standard library of theoretical stellar spectra intended for
    multiple synthetic photometry applications including spectral
    evolutionary synthesis is presented. The grid includes M dwarf model
    spectra, hence complementing the first library version established in
    Paper I (Lejeune et al., 1997, Cat. <J/A+AS/125/229>). It covers the
    following wide ranges of fundamental parameters:

      Teff: 50,000 to 2000K,
      logg:    5.5 to -1.02, and
    [Fe/H]:   +1.0 to -5.0.

    A correction procedure is also applied to the theoretical spectra in
    order to provide color-calibrated flux distributions over a large
    domain of effective temperatures. Empirical Teff-color calibrations
    are constructed between 11500K and 2000K, and semi-empirical
    calibrations for non-solar abundances ([Fe/H]=-3.5 to +1.0) are
    established. Model colors and bolometric corrections for both the
    original and the corrected spectra, synthesized in the
    UBV(RI)c(JHKLL'M) system, are given for the full range of stellar
    parameters.

    Synthetic colors:
    ----------------
    Synthetic UBV(RI)c(JHKLM) colors have been computed from both the
    original and the corrected model flux distributions presented in
    Paper I (1997A&AS..125..229L; see catalog <J/A+AS/125/229>), as
    the files lcb98ori.dat and lcb98cor.dat respectively; the results
    are also presented in individual files lcb98xxx.ori and lcb98xxx.cor,
    where xxx designates the metallicity (ex: 'm15' --> [Fe/H]=-1.5).
    For each file, we give synthetic colors computed from energy-weighted
    and photon-weighted stellar fluxes.

    Semi-empirical calibrations:
    ---------------------------

    Empirical ([Fe/H]=0.0) and semi-empirical (-3.5<=[Fe/H]<=+1.0)
    Teff-colors (UBVRIJHKLM) calibrations are given in Tables 1 to 10.

See also:
   VI/39 : Model Atmospheres (Kurucz, 1979)
   VI/78 : Theoretical Stellar Flux Spectra for F- to K-type Stars (Buser+ 1992)
   J/A+AS/105/311 : M giants spectra and photometry (Fluks+, 1994)
   J/A+AS/125/229 : A standard stellar library (Lejeune+ 1997)

File Summary:
--------------------------------------------------------------------------------
  FileName   Lrecl  Records   Explanations
--------------------------------------------------------------------------------
ReadMe         80      .   This file
table1.dat    100     47   [Fe/H]=0.0 empirical UBVRIJHKLM colors and
                                theoretical bolometric corrections
table2.dat    100     47   [Fe/H]=-3.5 empirical UBVRIJHKLM colors and
                                theoretical bolometric corrections
table3.dat    100     47   [Fe/H]=-3.0 empirical UBVRIJHKLM colors and
                                theoretical bolometric corrections
table4.dat    100     42   [Fe/H]=-2.5 empirical UBVRIJHKLM colors and
                                theoretical bolometric corrections
table5.dat    100     42   [Fe/H]=-2.0 empirical UBVRIJHKLM colors and
                                theoretical bolometric corrections
table6.dat    100     42   [Fe/H]=-1.5 empirical UBVRIJHKLM colors and
                                theoretical bolometric corrections
table7.dat    100     42   [Fe/H]=-1.0 empirical UBVRIJHKLM colors and
                                theoretical bolometric corrections
table8.dat    100     42   [Fe/H]=-0.5 empirical UBVRIJHKLM colors and
                                theoretical bolometric corrections
table9.dat    100     47   [Fe/H]=+0.5 empirical UBVRIJHKLM colors and
                                theoretical bolometric corrections
table10.dat   100     35   [Fe/H]=+1.0 empirical UBVRIJHKLM colors and
                                theoretical bolometric corrections
tables.dat    100    433   All empirical UBVRIJHKLM colors and theoretical
                                bolometric corrections (tables 1-10)
lcb98ori.dat  266   8316   All files lcb98*.ori: Synthetic colors from
                                original Stellar Energy Distributions (SEDS)
lcb98cor.dat  266   8315   All files lcb98*.cor: Synthetic colors from
                                corrected Stellar Energy Distributions (SEDS)
lcb98m01.ori  266    453   Synthetic colors from original SEDS at [Fe/H] = -0.1
lcb98m02.ori  266    453   Synthetic colors from original SEDS at [Fe/H] = -0.2
lcb98m03.ori  266    453   Synthetic colors from original SEDS at [Fe/H] = -0.3
lcb98m05.ori  266    456   Synthetic colors from original SEDS at [Fe/H] = -0.5
lcb98m10.ori  266    454   Synthetic colors from original SEDS at [Fe/H] = -1.0
lcb98m15.ori  266    456   Synthetic colors from original SEDS at [Fe/H] = -1.5
lcb98m20.ori  266    429   Synthetic colors from original SEDS at [Fe/H] = -2.0
lcb98m25.ori  266    427   Synthetic colors from original SEDS at [Fe/H] = -2.5
lcb98m30.ori  266    436   Synthetic colors from original SEDS at [Fe/H] = -3.0
lcb98m35.ori  266    430   Synthetic colors from original SEDS at [Fe/H] = -3.5
lcb98m40.ori  266    423   Synthetic colors from original SEDS at [Fe/H] = -4.0
lcb98m45.ori  266    393   Synthetic colors from original SEDS at [Fe/H] = -4.5
lcb98m50.ori  266    387   Synthetic colors from original SEDS at [Fe/H] = -5.0
lcb98p00.ori  266    467   Synthetic colors from original SEDS at [Fe/H] =  0.0
lcb98p01.ori  266    461   Synthetic colors from original SEDS at [Fe/H] =  0.1
lcb98p02.ori  266    457   Synthetic colors from original SEDS at [Fe/H] =  0.2
lcb98p03.ori  266    454   Synthetic colors from original SEDS at [Fe/H] =  0.3
lcb98p05.ori  266    452   Synthetic colors from original SEDS at [Fe/H] =  0.5
lcb98p10.ori  266    375   Synthetic colors from original SEDS at [Fe/H] =  1.0
lcb98m01.cor  266    453   Synthetic colors from corrected SEDS at [Fe/H] = -0.1
lcb98m02.cor  266    453   Synthetic colors from corrected SEDS at [Fe/H] = -0.2
lcb98m03.cor  266    453   Synthetic colors from corrected SEDS at [Fe/H] = -0.3
lcb98m05.cor  266    456   Synthetic colors from corrected SEDS at [Fe/H] = -0.5
lcb98m10.cor  266    454   Synthetic colors from corrected SEDS at [Fe/H] = -1.0
lcb98m15.cor  266    456   Synthetic colors from corrected SEDS at [Fe/H] = -1.5
lcb98m20.cor  266    429   Synthetic colors from corrected SEDS at [Fe/H] = -2.0
lcb98m25.cor  266    426   Synthetic colors from corrected SEDS at [Fe/H] = -2.5
lcb98m30.cor  266    436   Synthetic colors from corrected SEDS at [Fe/H] = -3.0
lcb98m35.cor  266    430   Synthetic colors from corrected SEDS at [Fe/H] = -3.5
lcb98m40.cor  266    423   Synthetic colors from corrected SEDS at [Fe/H] = -4.0
lcb98m45.cor  266    393   Synthetic colors from corrected SEDS at [Fe/H] = -4.5
lcb98m50.cor  266    387   Synthetic colors from corrected SEDS at [Fe/H] = -5.0
lcb98p00.cor  266    467   Synthetic colors from corrected SEDS at [Fe/H] =  0.0
lcb98p01.cor  266    461   Synthetic colors from corrected SEDS at [Fe/H] =  0.1
lcb98p02.cor  266    457   Synthetic colors from corrected SEDS at [Fe/H] =  0.2
lcb98p03.cor  266    454   Synthetic colors from corrected SEDS at [Fe/H] =  0.3
lcb98p05.cor  266    452   Synthetic colors from corrected SEDS at [Fe/H] =  0.5
lcb98p10.cor  266    375   Synthetic colors from corrected SEDS at [Fe/H] =  1.0
--------------------------------------------------------------------------------

Byte-by-byte Description of file: table*.dat
--------------------------------------------------------------------------------
   Bytes Format Units   Label     Explanations
--------------------------------------------------------------------------------
       1  A1    ---     Type      [dg] d: dwarf; g: giant
   3-  8  F6.0  K       Teff      Effective temperature
   9- 15  F7.3  mag     U-B       (U-B) colour
  17- 23  F7.3  mag     B-V       (B-V) colour
  25- 31  F7.3  mag     V-I       (V-I) colour
  33- 39  F7.3  mag     V-K       (V-K) colour
  41- 47  F7.3  mag     R-I       (R-I) colour
  49- 55  F7.3  mag     J-H       (J-H) colour
  57- 63  F7.3  mag     H-K       (H-K) colour
  65- 71  F7.3  mag     J-K       (J-K) colour
  73- 79  F7.3  mag     K-L       (K-L) colour
  81- 87  F7.3  mag     BCVmag    Bolometric correction in V (energy-weighted)
  89- 95  F7.3  ---     logg      Surface gravity
  97-100  F4.1 [Sun]    [Fe/H]    Metallicity
--------------------------------------------------------------------------------

Byte-by-byte Description of file: lcb98*.cor lcb98*.ori lcb*.dat
--------------------------------------------------------------------------------
   Bytes Format Units   Label     Explanations
--------------------------------------------------------------------------------
   1-  6  F6.0  K       Teff      Effective temperature
   7- 12  F6.2  ---     logg      Surface gravity
  13- 18  F6.2  [Sun]   [Fe/H]    Metallicity
  19- 26  F8.3  mag     BOL       Bolometric magnitude
  27- 34  F8.3  mag     BCVmag    Bolometric correction in V (energy-weighted)
  36- 42  F7.3  mag     U-B       (U-B) colour (energy-weighted)
  44- 50  F7.3  mag     B-V       (B-V) colour (energy-weighted)
  52- 58  F7.3  mag     V-R       (V-R) colour (energy-weighted)
  60- 66  F7.3  mag     V-I       (V-I) colour (energy-weighted)
  68- 74  F7.3  mag     V-K       (V-K) colour (energy-weighted)
  76- 82  F7.3  mag     R-I       (R-I) colour (energy-weighted)
  84- 90  F7.3  mag     I-K       (I-K) colour (energy-weighted)

  92- 98  F7.3  mag     J-H       (J-H) colour (energy-weighted)
 100-106  F7.3  mag     H-K       (H-K) colour (energy-weighted)
 108-114  F7.3  mag     K-L       (K-L) colour (energy-weighted)
 116-122  F7.3  mag     J-K       (J-K) colour (energy-weighted)
 124-130  F7.3  mag     J-L       (J-L) colour (energy-weighted)
 132-138  F7.3  mag     J-L'      (J-L')colour (energy-weighted)
 140-146  F7.3  mag     K-M       (K-M) colour (energy-weighted)
 147-154  F8.3  mag     BCVmagp   Bolometric correction in V (photon-weighted)
 156-162  F7.3  mag     U-Bp      (U-B) colour (photon-weighted)
 164-170  F7.3  mag     B-Vp      (B-V) colour (photon-weighted)
 172-178  F7.3  mag     V-Rp      (V-R) colour (photon-weighted)
 180-186  F7.3  mag     V-Ip      (V-I) colour (photon-weighted)
 188-194  F7.3  mag     V-Kp      (V-K) colour (photon-weighted)
 196-202  F7.3  mag     R-Ip      (R-I) colour (photon-weighted)
 204-210  F7.3  mag     I-Kp      (I-K) colour (photon-weighted)
 212-218  F7.3  mag     J-Hp      (J-H) colour (photon-weighted)
 220-226  F7.3  mag     H-Kp      (H-K) colour (photon-weighted)
 228-234  F7.3  mag     K-Lp      (K-L) colour (photon-weighted)
 236-242  F7.3  mag     J-Kp      (J-K) colour (photon-weighted)
 244-250  F7.3  mag     J-Lp      (J-L) colour (photon-weighted)
 252-258  F7.3  mag     J-L'p     (J-L')colour (photon-weighted)
 260-266  F7.3  mag     K-Mp      (K-M) colour (photon-weighted)
--------------------------------------------------------------------------------

Acknowledgements: Thibault Lejeune <lejeune@astro.unibas.ch>
================================================================================
(End)                                          Patricia Bauer [CDS]  18-Jun-1998
