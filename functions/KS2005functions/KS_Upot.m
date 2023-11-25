function U = KS_Upot(x, y, z, a, c, p, r, M)
% Calculate magnetic potential terms from model coefficients for magnetotail field contribution.
%
% Parameters
% ----------
% x : double, 1xN
%   x component of evaluation point in the System III frame.
% y : double, 1xN
%   y component of evaluation point in the System III frame.
% z : double, 1xN
%   z component of evaluation point in the System III frame.
% a : double, 8x64
%   Magnetotail parameterization coefficients.
% c : double, 8x64
%   Magnetotail parameterization coefficients.
% p : double, 8x16
%   Magnetotail parameterization coefficients.
% r : double, 8x16
%   Magnetotail parameterization coefficients.
% M : int
%   Dimension (MxM) of magnetotail coefficients ``a`` and ``c``, typically ``Mode + 1``.
%
% Returns
% -------
% outputName : type, dims
%   Description.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Initialize output
    U = zeros(size(x));

    M2 = 2*M;
    ix = 1;
    
    % Calculate a/p terms and add into U
    for i=1:M
        for k=(M+1):M2
            T = exp(sqrt(1/p(i)^2 + 1/p(k)^2) * x) .* cos(y/p(i)) .* sin(z/p(k));
            U = U + a(i, k-M) * T;
            ix = ix + 1;
        end
    end
    
    % Calculate c/r terms and sum with a/p terms
    for i=1:M
        for k=(M+1):M2
            T = exp(sqrt(1/r(i)^2 + 1/r(k)^2) * x) .* sin(y/r(i)) .* sin(z/r(k));
            U = U + c(i, k-M) * T;
            ix = ix + 1;
        end
    end
    
end
