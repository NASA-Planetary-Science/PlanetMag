Implemented frames and coordinates
==================================

Numerous frames are implemented in SPICE through built-in options and those we have defined in ``custom_frames_v01.tf``. We describe each implemented frame and its associate coordinate system here. More detail on the frames and where they have appeared in past studies can be found in the supplemental text of the PlanetMag publication (DOI TBD).

All frames use the body center of mass for the origin. Any axis direction not listed can be inferred by that needed to complete a right-handed coordinate system. All frames in Planet\Mag and all those available in SPICE are right-handed. All frames implemented in Planet\Mag use Cartesian basis vectors.

IAU frames
++++++++++

Defined in reports of the IAU Working Group on Cartographic Coordinates and Rotational Elements.

* :math:`\hat{z}` is along the angular momentum vector
* :math:`\hat{x}` is approximately toward the planet, with the :math:`xz` plane containing a specific surface feature identified by the IAU
* Available in SPICE as ``IAU_BODY`` for any body for which ephemeris data has been loaded into the kernel pool

System III frames
+++++++++++++++++

These frames rotate with the planet according to the period inferred from regular oscillations in the magnetic field, or physical behavior tied to the field such as synchrotron emission.

* :math:`\hat{z}` is along the planet's angular momentum vector
* :math:`\hat{x}` is in an arbitrarily defined direction at an arbitrary moment, because the giant planets generally lack stable surface features
* Identical to ``IAU_JUPITER`` and ``IAU_SATURN`` for those planets
* Implemented as ``ULS`` and ``NLS`` for the Uranus Longitude System and Neptune Longitude System, respectively

SPRH frames
+++++++++++

This is a designation sometimes used in NASA Planetary Data System (PDS) archives for data aligned to spherical coordinate basis vectors. The frame is identical to System III coordinates for those bodies for which data are expressed in this format. However, because spherical coordinate axes are not aligned to Cartesian axes, the data must be converted for comparison to SPICE frames.

Planet--Sun--Orbit frames
+++++++++++++++++++++++++

* :math:`\hat{x}` is toward the Sun
* :math:`\hat{y}` is along the component normal to :math:`\hat{x}` of the velocity vector of the Sun as seen from the body
* Implemented as ``JSO``, ``KSO``, ``USO``, ``NSO`` for Jupiter, Saturn, Uranus, and Neptune respectively

Planet--Sun--Magnetic frames
++++++++++++++++++++++++++++

* :math:`\hat{x}` is toward the Sun
* :math:`\hat{y}` is along :math:`\hat{\mathbf{M}}\times\hat{x}`, where :math:`\mathbf{M}` is the magnetic dipole moment vector. Requires a choice of magnetic field model in order to infer the dipole moment direction.
* For use in particular models that require this frame. The selected multipole models are:

   * ``JSM``: O4 (Acuna and Ness, 1976) https://doi.org/10.1029/JA081i016p02917
   * ``KSM``: Cassini 11 (Dougherty et al., 2018) https://doi.org/10.1126/science.aat5434
   * ``USM``: OTD (Ness et al., 1986) https://doi.org/10.1126/science.233.4759.85
   * ``NSM``: O8 (Connerney et al., 1992) https://doi.org/10.1016/0273-1177(92)90394-D

Planet--Dipole--Solar--Zenith frames
++++++++++++++++++++++++++++++++++++

* :math:`\hat{z}` is toward the Sun
* :math:`\hat{y}` is along :math:`\hat{\mathbf{M}}\times\hat{z}`, where :math:`\mathbf{M}` is the same dipole moment vector as used in Planet--Sun--Magnetic frames
* Implemented as ``JDSZ``, ``KDSZ``, ``UDSZ``, ``NDSZ`` for Jupiter, Saturn, Uranus, and Neptune respectively

Solar--Magnetic--Planet frames
++++++++++++++++++++++++++++++

* :math:`\hat{z}` is along :math:`\hat{\mathbf{M}}` as defined in Planet--Sun--Magnetic frames
* :math:`\phi=0` at the sub-solar point, i.e. the :math:`x` axis is along the direction perpendicular to :math:`\hat{z}` along the noon local time longitude line
* The :math:`xy` plane defines the magnetic equatorial plane
* Implemented as ``SMJ``, ``SMK``, ``SMU``, ``SMN`` for Jupiter, Saturn, Uranus, and Neptune respectively
