# PlanetMag
<img src="misc/PlanetMag_logo.png" width=150 title="PlanetMag logo"/>

Matlab software for evaluation of planetary magnetic field models and the oscillations applied to moons. PlanetMag uses the [SPICE toolkit](https://naif.jpl.nasa.gov/naif/toolkit.html) and published magnetospheric models derived from spacecraft data. All timing and phases are referenced to the J2000 epoch, 2000-Jan-01 12:00:00 TDB, which is equivalent to 11:58:55.816 am UTC on the same date.

The main repository is mirrored at https://github.com/NASA-Planetary-Science/PlanetMag; any pull requests should be submitted to https://github.com/coreyjcochrane/PlanetMag.

## Installation
1. PlanetMag requires Mice — [the Matlab implementation of the SPICE toolkit](https://naif.jpl.nasa.gov/naif/toolkit_MATLAB.html) — and specific SPICE kernels. Mice functions must be located on the Matlab path. The specific SPICE kernels used are all listed in `LoadSpice.m`, and the expected path to find them is in the same file (`./spice/` by default). Smaller kernels are packaged with PlanetMag, but larger generic kernels must be downloaded from the [generic kernels for satellites](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/) page of the NAIF website. Earth data requires the latest Earth PCK file available at https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/earth_latest_high_prec.bpc.
1. Some mission data (especially from ESA) uses Common Data Format (.cdf) files. Reading these files requires the CDF toolkit, available at https://spdf.gsfc.nasa.gov/pub/software/cdf/dist/latest/. Matlab ships with a version of the toolkit that works for most files, but can lock up reading large files (>15 MB).
    1. An example script for converting CDF files to ASCII tables is in ConvertSwarmCDF.m.

## Running the software
1. Run `PlanetMag.m` from the directory it is in to test for full functionality and a demonstration of evaluation of excitation moments (amplitude and phase of oscillations) applied to Europa by Jupiter.
1. Comparison between magnetic models (currently only implemented for Jupiter) requires downloading spacecraft data from PDS into MAG/Mission, where "Mission" is the name of the spacecraft from which data will be compared. .tab files from PDS should keep the same filenames. 
    1. Swarm data is used for Earth field model validation. Data files are available at https://swarm-diss.eo.esa.int/#swarm/Level1b/Latest_baselines/MAGx_LR 

## Acknowledging PlanetMag
We appreciate your interest in PlanetMag! Please consider alerting us to your work (corey.j.cochrane@jpl.nasa.gov). Suggested acknowledgement in publications: "Data used in this work were generated using the open source PlanetMag framework hosted on GitHub."
