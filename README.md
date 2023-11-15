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
1. Download mice and place it in the right spot
1. Download the following SPICE kernels to the ``spice`` directory:
    1. kernel 1
    1. kernel 2
1. Download spacecraft magnetic data for comparison to ``MAG`` directory

## Running the software
1. Run `PlanetMag.m` from the directory it is in to test for full functionality and a demonstration of evaluation of excitation moments (amplitude and phase of oscillations) applied to Europa by Jupiter.
1. Comparison between magnetic models requires downloading spacecraft data from PDS into MAG/Mission, where "Mission" is the name of the spacecraft from which data will be compared. .tab files from PDS should keep the same filenames. 
    1. Compare field models evaluated with `PlanetMag.m` using the `Comparison` scripts in the `comparison` directory.
    1. Swarm data is used for Earth field model validation. Data files are available at <https://swarm-diss.eo.esa.int/#swarm/Level1b/Latest_baselines/MAGx_LR>

## Acknowledging PlanetMag
We appreciate your interest in PlanetMag! Please consider alerting us to your work (corey.j.cochrane@jpl.nasa.gov). Suggested acknowledgement in publications: "Data used in this work were generated using the open-source _PlanetMag_ framework hosted on GitHub (<https://github.com/coreyjcochrane/PlanetMag>)."

## Copyright notice
Copyright 2023, by the California Institute of Technology. ALL RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any commercial use must be negotiated with the Office of Technology Transfer at the California Institute of Technology.

This software may be subject to U.S. export control laws. By accepting this software, the user agrees to comply with all applicable U.S. export laws and regulations. User has the responsibility to obtain export licenses, or other export authority as may be required before exporting such information to foreign countries or providing access to foreign persons.
