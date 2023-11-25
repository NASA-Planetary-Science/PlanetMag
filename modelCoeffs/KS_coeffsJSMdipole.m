function dipole = KS_coeffsJSMdipole
% Import parameters describing the orientation of the magnetic dipole moment in the JSM frame.
%
% Includes a 202 degree rotation about z and a 9.6 degree rotation about y. According to the
% following document, this functions is therefore using "outdated" O4 model parameters: 
% https://lasp.colorado.edu/home/mop/files/2015/02/CoOrd_systems12.pdf
%
% Returns
% -------
% dipole : double, 3x3
%   Rotation matrix used in converting between the System III frame and the Jupiter--Sun--magnetic
%   frame.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    dipole = [-0.9141996,   0.36936062,  -0.16676875; ...
              -0.3746066,  -0.92718385,            0; ...
              -0.154625290, 0.0624726744, 0.98599604];
end
