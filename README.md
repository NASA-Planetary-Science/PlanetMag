# PlanetMag
![PlanetMag logo](misc/PlanetMag_logoDocs.png)

Matlab software for evaluation of planetary magnetic field models and the oscillations applied to moons. The primary purposes of this software are to offer capabilities for:
* Evaluation of planetary magnetic fields at arbitrary locations and times around the supported planets based on models described in peer-reviewed publications based on spacecraft data. The supported planets are Earth, Jupiter, Saturn, Uranus, and Neptune.
* Determination of frequency spectra of magnetic oscillations in the frame of reference of satellites orbiting each supported planet, including recording complex amplitudes in each vector component in standard coordinate systems referenced to the J2000.0 epoch.
* Statistical comparison of planetary field model results against spacecraft data, for validation purposes.

PlanetMag uses the [SPICE toolkit](https://naif.jpl.nasa.gov/naif/toolkit.html) and published magnetospheric models derived from spacecraft data. All timing and phases are referenced to the J2000 epoch, 2000-Jan-01 12:00:00 TDB, which is equivalent to 11:58:55.816 am UTC on the same date.

The main repository is mirrored at <https://github.com/NASA-Planetary-Science/PlanetMag>; any pull requests should be submitted to the primary repository at <https://github.com/coreyjcochrane/PlanetMag>. Read the software documentation at <https://coreyjcochrane.github.io/PlanetMag>.

This scientific software is the brain child of Corey J. Cochrane and Erik Sturm. The PlanetMag framework was developed by Corey J. Cochrane and Marshall J. Styczinski. The software is currently maintained by Marshall J. Styczinski.

Questions about PlanetMag? Please contact Corey J. Cochrane at corey.j.cochrane@jpl.nasa.gov. If you run into problems or need help with using the software, please contact maintainer Marshall J. Styczinski at itsmoosh@gmail.com.

## Installation
PlanetMag requires Mice — [the Matlab implementation of the SPICE toolkit](https://naif.jpl.nasa.gov/naif/toolkit_MATLAB.html) — and specific SPICE kernels. Mice functions must be located on the Matlab path. The specific SPICE kernels used are all listed in `LoadSpice.m`, and the expected path to find them is in the same file (`./spice/` by default). Smaller kernels are packaged with PlanetMag, but larger generic kernels must be downloaded from the [generic kernels for satellites](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/) page of the NAIF website. Earth data requires the latest Earth PCK file available at <https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/earth_latest_high_prec.bpc>.

Some mission data (especially from ESA) uses Common Data Format (.cdf) files. Reading these files requires the CDF toolkit, available at <https://spdf.gsfc.nasa.gov/pub/software/cdf/dist/latest/>. Matlab ships with a version of the toolkit that works for most files, but can lock up reading large files (>15 MB). An example script for converting CDF files to ASCII tables is in ConvertSwarmCDF.m.

To fully install PlanetMag and access all features:
1. Download or clone Matlab repo
1. Download Mice and place it in the right spot
1. Download the following SPICE kernels to the ``spice`` directory:
   1. ``de430.bsp`` --- from <https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/>
   1. ``jup365.bsp`` --- from <https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/>
   1. ``sat441.bsp`` --- from <https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/>
   1. ``ura111.bsp`` --- from <https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/>
   1. ``nep097.bsp`` --- from <https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/>
   1. ``earth_latest_high_prec.bpc`` --- from <https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/>
   1. For Galileo, all in the ``spice/Galileo`` directory and all from <https://naif.jpl.nasa.gov/pub/naif/GLL/kernels/spk/>:
         1. ``s980326a.bsp``
         1. ``s000131a.bsp``
         1. ``s030916a.bsp``
   1. For Juno, all from <https://naif.jpl.nasa.gov/pub/naif/JUNO/kernels/spk/>:
         1. In the ``spice`` directory:
            1.  ``jup380s.bsp``
         1. In the ``spice/Juno`` directory:
            1.  ``juno_rec_orbit.bsp``
   1. For Voyager 1, all in the ``spice/Voyager 1`` directory and all from <https://naif.jpl.nasa.gov/pub/naif/VOYAGER/kernels/spk/>:
         1. ``vgr1_jup230.bsp``
         1. ``vgr1_sat337.bsp``
         1. ``s030916a.bsp``
   1. For Voyager 2, all in the ``spice/Voyager 2`` directory and all from <https://naif.jpl.nasa.gov/pub/naif/VOYAGER/kernels/spk/>:
         1. ``vgr2_jup230.bsp``
         1. ``vgr2_sat337.bsp``
         1. ``vgr2.ura111.bsp``
         1. ``vgr2_nep097.bsp``
1. Download spacecraft magnetic data for comparison to the ``MAG`` directory

## Running the software
1. Run `PlanetMag.m` from the directory it is in to test for full functionality and a demonstration of evaluation of excitation moments (amplitude and phase of oscillations) applied to Europa by Jupiter.
1. Comparison between magnetic models requires downloading spacecraft data from PDS into MAG/Mission, where "Mission" is the name of the spacecraft from which data will be compared. .tab files from PDS should keep the same filenames. 
    1. Compare field models evaluated with `PlanetMag.m` using the `Comparison` scripts in the `comparison` directory.
    1. Swarm data is used for Earth field model validation. Data files are available at <https://swarm-diss.eo.esa.int/#swarm/Level1b/Latest_baselines/MAGx_LR>

## Acknowledging PlanetMag
We appreciate your interest in PlanetMag! Please consider alerting us to your work (corey.j.cochrane@jpl.nasa.gov). Suggested acknowledgement in publications: "Data used in this work were generated using the open-source _PlanetMag_ framework hosted on GitHub (<https://github.com/coreyjcochrane/PlanetMag>)."

## Copyright notice
Copyright 2023, by the California Institute of Technology. ALL RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any commercial use must be negotiated with the Office of Technology Transfer at the California Institute of Technology.

## Disclaimer
This software may be subject to U.S. export control laws. By accepting this software, the user agrees to comply with all applicable U.S. export laws and regulations. User has the responsibility to obtain export licenses, or other export authority as may be required before exporting such information to foreign countries or providing access to foreign persons. This software is provided "as is" without any warranty, express or implied, and the recipient waives any claims against, and indemnifies and holds harmless, NASA and its contractors and subcontractors.
