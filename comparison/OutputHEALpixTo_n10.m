% This script prints an evaluation of Legendre polynomials and their
% derivatives on a HEALpix map, for validation purposes. The HEALpix
% map must be supplied in a text file at the specified location. All
% harmonics are mapped, from n=1 to n=10.

% HEALpix file information
hpDir = fullfile('~', 'MoonMag', 'outData');
hpLocsFile = fullfile(hpDir, 'healpix_locs.txt');

% Output file information
outDir = 'out';
fnameBase = fullfile(outDir, 'PM_Legendre_');

% Load HEALpix pixel locations from file
hpLocs = readtable(hpLocsFile, 'NumHeaderLines', 2);
theta_rad = hpLocs.Var1';
phi_rad = hpLocs.Var2';
r_RP = ones(size(theta_rad));

% Assign common strings
pureHarmFile = fullfile('modelCoeffs', 'coeffsPureHarmonic.csv');
pureHarmHeader = 'Contains ''g'' and ''h'' Schmidt semi-normalized spherical harmonic coefficients for evaluation of a pure harmonic. For testing and validation purposes.';
pureHarmCols = sprintf('%12s,%12s,%12s,%12s', 'n', 'm', 'gnm', 'hnm');
outHarmFileBase = fullfile('out', 'pureHarmMap_');

% Common value for g/h
ghVal = 1.0;

% Loop over n and m for g and h
for n = 1:10
    for m = 0:n
        WritePureHarmFile(pureHarmFile, pureHarmHeader, pureHarmCols, outHarmFileBase, ...
            n, m, ghVal, 0.0, r_RP, theta_rad, phi_rad, 'g');
        if m ~= 0
            WritePureHarmFile(pureHarmFile, pureHarmHeader, pureHarmCols, outHarmFileBase, ...
                n, m, 0.0, ghVal, r_RP, theta_rad, phi_rad, 'h');
        end
    end
end


function WritePureHarmFile(fpath, header, colNames, outFileBase, ...
    n, m, gnm, hnm, r_RP, theta_rad, phi_rad, outFlag)
    % Print a formatted file containing the specified pure harmonic.
    coeffsLine = sprintf('%12d,%12d,%12.8f,%12.8f', n, m, gnm, hnm);
    dlmwrite(fpath, header, 'delimiter', '');
    dlmwrite(fpath, colNames, 'delimiter','', '-append');
    dlmwrite(fpath, coeffsLine, 'delimiter','', 'precision',8, '-append');

    [Bvec, ~, ~] = MagFldParent('PureHarmonic', r_RP, theta_rad, phi_rad, 'SphericalHarmonic', 'None', 0, 0);
    save([outFileBase sprintf('n%02dm%02d%s.mat', n, m, outFlag)], 'r_RP', 'theta_rad', 'phi_rad', 'Bvec');
end