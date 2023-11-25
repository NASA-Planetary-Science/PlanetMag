function ptInsideMpause = KS_CheckIfInsideMP(xJSM, yJSM, zJSM)
% Return logical array indicating whether each point is inside the unmapped magnetopause.
%
% Parameters
% ----------
% xJSM : double, 1xN
%   x coordinate of evaluation point in Jupiter--Sun--magnetic (JSM) frame in planetary radii. See
%   KS_BJSMtoBxyz for a definition of the JSM frame.
% yJSM : double, 1xN
%   y coordinate of evaluation point in JSM frame in planetary radii.
% zJSM : double, 1xN
%   z coordinate of evaluation point in JSM frame in planetary radii.
%
% Returns
% -------
% ptInsideMpause : bool, 1xN
%   Whether each evaluation point is inside (true) or outside (false) of the magnetopause.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % MJS note: In K. Khurana's code, this function is split into day- and night-side calculations,
    % but everything being done for each has identical results. There are even typos copy-pasted
    % between the two blocks (Fortran is not case-sensitive so the typos are harmless here). I
    % don't see how or why there is any separate handling there.

    % Initialize
    ptInsideMpause = ones(size(xJSM));
    % Dynamic pressure
    DP = 0.04;
    % Symmetric magnetopause parameters
    A =   13779.154626782700000;
	B =    -130.046991164773700;
	C = -2.216972473764874D-001;
	D =  0.000000000000000D+000;
	E = -8.453162863101464D-001;
	F =  0.000000000000000D+000;
    
    phi = atan2(zJSM, yJSM);
    rho_in = sqrt(zJSM.^2 + yJSM.^2);
    
    % Calculate numbers needed to determine magnetopause distance
    aa = E.*cos(phi).^2 - sin(phi).^2 + F.*cos(phi).*sin(phi);
    bb = D.*cos(phi);
    cc = A + B.*xJSM + C.*xJSM.^2;
    term = bb.^2 - 4*aa.*cc;
    rho1 = (-bb + sqrt(term)) ./ (2*aa);
    rho2 = (-bb - sqrt(term)) ./ (2*aa);
    rho = zeros(size(aa));
    rho(rho2 > rho1) = rho2(rho2 > rho1);
    rho(~(rho2 > rho1)) = rho1(~(rho2 > rho1));
    
    % Set result for bounding terms
    ptInsideMpause(term < 0) = 0;
    ptInsideMpause(rho < 0) = 0;
    ptInsideMpause(rho_in > rho) = 0;
end
