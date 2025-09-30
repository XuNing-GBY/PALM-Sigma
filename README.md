# User’s Guide

Xu Ning: xuning.0912@gmail.com


## I. Introduction

Code Name: PALM-Sigma

Code Version: v1.0

Description: PALM-Sigma consists of a series of user-defined modules and modified PALM modules, developed for wave-phase-resolved LES. PALM-Sigma should be used as USER\_CODE based on the framework of PALM.

License: GNU General Public License v3.0

If you use PALM-Sigma in your paper, please cite the following paper (detail will be updated in PALM-Sigma's [GitHub's repository](https://github.com/XuNing-GBY/PALM-Sigma.git) when the paper is published):

>  Xu Ning, Mostafa Bakhoday-Paskyabi, 2025: Implementation of a sigma coordinate system in PALM-Sigma v1.0 (based on PALM v21.10) for LES study of the marine atmospheric boundary layer. *Submitted to Geoscientific Model Development*.

and the software as indicated in PALM-Sigma's Zenodo record.


## II. Software Requirements

PALM-Sigma v1.0 is built on [PALM v21.10](https://gitlab.palm-model.org/releases/palm\_model\_system). For instructions on installing and compiling PALM, please refer to the documentation available on the [PALM's website](https://palm.muk.uni-hannover.de/trac/wiki). Running PALM-Sigma based on other PALM versions may not work properly.


## III. Compilation & Running

All source files of PALM-Sigma should be included in the user-specific interface folder, i.e., USER\_CODE, which should be placed under the case folder together with the INPUT folder.

Compliation will be done automatically when a PALM job is submitted.


## IV. Configuration

PALM-Sigma uses the same INPUT files as in PALM, but with an extra NAMELIST "wind\_wave\_inta\_parameters".

The descriptions of five parameters are as following:

wwinta\_mode: CHARACTER type, the mode for wind-wave interaction, should always be 'sfc'.

wave\_type: CHARACTER type, the type of wave, possible options are 'regular' (regular wave that can be defined by a 2D sinusoidal function), and 'steady' (fixed bumps without moving).

wave\_amplitude: REAL type, wave amplitude, unit meter.

wave\_number: REAL type, wave number.

wave\_direction: REAL type, wave direction, unit degree. The wave direction is defined as the direction from which the waves propagate. North is set to 0.0°, and the angle increases clockwise. Thus, waves coming from the west correspond to 270.0°.

Notice that a too large wave height may cause a crash, a wave steepness (i.e., wave amplitude * wave number) less than 0.2 would be recommended.

