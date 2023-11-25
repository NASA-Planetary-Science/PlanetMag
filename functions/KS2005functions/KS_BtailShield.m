function [BxJSO, ByJSO, BzJSO] = KS_BtailShield(M, Mode, xJSO, yJSO, zJSO)
% Get magnetotail contribution to shielded dipole field.
%
% See KS_S3CtoJSO for a definition of the Jupiter--Sun--Orbital (JSO) frame.
%
% Parameters
% ----------
% M : int
%   Dimension (MxM) of magnetotail coefficients ``a`` and ``c``, typically ``Mode + 1``.
% Mode : int
%   Parameter selection for magnetotail model. Passing 7 or greater will include all coefficients.
% xJSO : double, 1xN
%   x coordinate of evaluation points in JSO frame.
% yJSO : double, 1xN
%   y coordinate of evaluation points in JSO frame.
% zJSO : double, 1xN
%   z coordinate of evaluation points in JSO frame.
%
% Returns
% -------
% BxJSO, ByJSO, BzJSO : double, 1xN
%   Magnetic field contribution from magnetotail in nT in JSO frame aligned with cartesian axes.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Initialize output
    [BxJSO, ByJSO, BzJSO] = deal(zeros(size(xJSO)));

    % Get loop counter from Mode input
    if Mode < 7
        Modes = Mode:Mode;
    else
        Modes = 1:6;
    end

    % Coefficient tables are scaled against distances in units of 100 RJ -- scale input coordinates
    % to match
    x = xJSO/100;
    y = yJSO/100;
    z = zJSO/100;
    
    xp = x + 0.001;
	xm = x - 0.001;
	yp = y + 0.001;
	ym = y - 0.001;
	zp = z + 0.001;
	zm = z - 0.001;
        
    % Fetch coefficients
    CS = KS_coeffsCsheet();
    [a, c, p, r] = KS_coeffsMtail();
    [ain, cin] = deal(zeros(M, M));
    
    % Collect a and c values for the specified mode(s)
    for thisMode=Modes
        icount = 0;
        for ii=1:M
            for kk=1:M
                icount = icount+1;
                ain(ii,kk) = a(thisMode,icount);
                cin(ii,kk) = c(thisMode,icount);
            end
        end
        
        % Calculate magnetic potentials and take the gradient to get Bvec
        Up = KS_Upot(xp, y, z, ain, cin, p(thisMode,:), r(thisMode,:), M);
        Um = KS_Upot(xm, y, z, ain, cin, p(thisMode,:), r(thisMode,:), M);
        Bx = -(Up - Um) / 0.002;
        
        Up = KS_Upot(x, yp, z, ain, cin, p(thisMode,:), r(thisMode,:), M);
        Um = KS_Upot(x, ym, z, ain, cin, p(thisMode,:), r(thisMode,:), M);
        By = -(Up - Um) / 0.002;
        
        Up = KS_Upot(x, y, zp, ain, cin, p(thisMode,:), r(thisMode,:), M);
        Um = KS_Upot(x, y, zm, ain, cin, p(thisMode,:), r(thisMode,:), M);
        Bz = -(Up - Um) / 0.002;
        
        Bx = Bx * CS(thisMode);
        By = By * CS(thisMode);
        Bz = Bz * CS(thisMode);
        
        BxJSO = BxJSO + Bx;
        ByJSO = ByJSO + By;
        BzJSO = BzJSO + Bz;
    end
    
    % Negate for output
    BxJSO = -BxJSO;
    ByJSO = -ByJSO;
    BzJSO = -BzJSO;
end
