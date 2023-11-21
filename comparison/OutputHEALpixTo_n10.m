function OutputHEALpixTo_n10(hpDir, hpFname, outDir, outFbase)
% Evaluates intermediate functions for each n,m on a HEALpix map and prints to disk.
%
% Prints an evaluation of Legendre polynomials and their derivatives on a HEALpix map, for
% validation purposes. The HEALpix map must be supplied in a text file at the specified location.
% All harmonics are mapped, from n=1 to n=10.
%
% Parameters
% ----------
% hpDir : char, 1xC, default=fullfile('~', 'MoonMag', 'outData')
%   Directory containing a file with pixel locations on a HEALpix grid where function outputs are
%   to be compared.
% hpFname : char, 1xD, default='healpix_locs.txt'
%   File name where HEALpix locations are stored.
% outDir : char, 1xE, default='out'
%   Directory where a file with function outputs will be written.
% outFbase : char, 1xF, default='pureHarmMap\_'
%   File name pattern for function outputs.
% ghVal : double, default=1.0
%   Common value for g and h in gauss used in MoonMag and P\lanetMag intermediate function
%   evaluation.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~exist('hpDir', 'var'); hpDir = fullfile('~', 'MoonMag', 'outData'); end
    if ~exist('hpFname', 'var'); hpFname = 'healpix_locs.txt'; end
    if ~exist('outDir', 'var'); outDir = 'out'; end
    if ~exist('outFbase', 'var'); outFbase = 'pureHarmMap_'; end
    if ~exist('ghVal', 'var'); ghVal = 1.0; end

    % HEALpix file information
    hpLocsFile = fullfile(hpDir, hpFname);

    % Load HEALpix pixel locations from file
    hpLocs = readtable(hpLocsFile, 'NumHeaderLines', 2);
    theta_rad = hpLocs.Var1';
    phi_rad = hpLocs.Var2';
    r_RP = ones(size(theta_rad));

    % Assign common strings
    pureHarmFile = fullfile('modelCoeffs', 'coeffsPureHarmonic.csv');
    pureHarmHeader = ['Contains ''g'' and ''h'' Schmidt semi-normalized spherical harmonic ' ...
        'coefficients for evaluation of a pure harmonic. For testing and validation purposes.'];
    pureHarmCols = sprintf('%12s,%12s,%12s,%12s', 'n', 'm', 'gnm', 'hnm');
    outHarmFileBase = fullfile(outDir, outFbase);

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

end
