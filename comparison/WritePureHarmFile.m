function WritePureHarmFile(fpath, header, colNames, outFileBase, n, m, gnm, hnm, r_RP, ...
    theta_rad, phi_rad, outFlag)
% Print a formatted file containing the specified pure harmonic.
%
% Parameters
% ----------
% fpath : char, 1xC
%   Full path to output file.
% header : char, 1xD
%   Description of the output file contents.
% colNames : char, 1xE
%   Formatted string describing output columns.
% outFileBase : char, 1xF
%   Pattern for output file names to which n and m are appended.
% n, m : int
%   Degree and order of output functions in Schmidt normalization.
% gnm, hnm : double
%   Value to use for Schmidt semi-normalized gauss coefficients in gauss.
% r_RP : double, 1xN
%   Radius in units of :math:`R_P` (planetary radii) at which to evaluate functions.
% theta_rad : double, 1xN
%   Colatitude in radians at which to evaluate functions.
% phi_rad : double, 1xN
%   East longitude in radians at which to evaluate functions.
% outFlag : char, 1xG
%   Strings to append to filenames to distinguish pure g from pure h for each n,m pair.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    coeffsLine = sprintf('%12d,%12d,%12.8f,%12.8f', n, m, gnm, hnm);
    dlmwrite(fpath, header, 'delimiter', '');
    dlmwrite(fpath, colNames, 'delimiter','', '-append');
    dlmwrite(fpath, coeffsLine, 'delimiter','', 'precision',8, '-append');

    [Bvec, ~, ~] = MagFldParent('PureHarmonic', r_RP, theta_rad, phi_rad, 'SphericalHarmonic', ...
                                'None', 0, 0);
    save([outFileBase sprintf('n%02dm%02d%s.mat', n, m, outFlag)], 'r_RP', 'theta_rad', ...
         'phi_rad', 'Bvec');
end
