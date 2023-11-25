function [g, h] = KS_coeffsVIP4(coeffPath, nHeadLines, AS_CODED)
% Import spherical harmonic coefficients for the VIP4 model.
%
% See GetGaussCoeffs for a more detailed description of the VIP4 model. This function includes the
% option to use less-precise values for the coefficients as originally implemented in the Khurana
% and Schwarzl (2005) model.
%
% Parameters
% ----------
% coeffPath : char, 1xE, default='modelCoeffs'
%   Directory where coefficients files are located.
% nHeadLines : int, default=2
%   Number of header lines in model coefficient files. This is typically ``2``: One for a
%   description of the file contents and one for column headers.
% AS_CODED : bool, default=0
%   Whether to match the original Fortran code (true) or with increased precision and corrected
%   parameters (false). See MagFldJupiterKS2005 for more details.
%
% Returns
% -------
% g, h : double, 4x5
%   Real spherical harmonic coefficients (Gauss coefficients) in Schmidt semi-normalized form for
%   the internal magnetic field of the planet in terms of surface field strength in gauss in the
%   VIP4 model.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~exist('coeffPath', 'var'); coeffPath = 'modelCoeffs'; end
    if ~exist('nHeadLines', 'var'); nHeadLines = 2; end

    if AS_CODED
        % This list matches the less-precise table used in K. Khurana's Jupiter field model code,
        % and all but the dipole moments match the table in the associated publication (the dipole
        % moments use more precise values).
        %        g0n,      g1n,      g2n,      g3n,     g4n
        g = [4.20543, -0.65920,        0,        0,       0;
            -0.05100, -0.61900,  0.49700,        0,       0;
            -0.01600, -0.52000,  0.24400, -0.17600,       0;
            -0.16800,  0.22200, -0.06100, -0.20200, 0.06600];
        %        h0n,      h1n,      h2n,      h3n,     h4n
        h =       [0,  0.24992,        0,        0,       0;
                   0, -0.36100,  0.05300,        0,       0;
                   0, -0.08800,  0.40800, -0.31600,       0;
                   0,  0.07600,  0.40400, -0.16600, 0.03900];
    else
        % More precise values as imported by GetGaussCoeffs
        g = dlmread(fullfile(coeffPath, 'coeffsJupiterVIP4g.csv'), ',', nHeadLines, 0);
        h = dlmread(fullfile(coeffPath, 'coeffsJupiterVIP4h.csv'), ',', nHeadLines, 0);
    end
end
