function [Bx, By, Bz] = MPboxHarmonic(xyzPSM_Rp, Psi, a, b, c, d, pqrs)
% Calculate the magnetic field contributuon at specific points based on a box harmonix model.
%
% Applies a box harnomic model based on Tsyganenko (2002) https://doi.org/10.1029/2001JA000219
% as formulated by Arridge and Eggington (2021) https://doi.org/10.1016/j.icarus.2021.114562.
%
% Parameters
% ----------
% xyzPSM_Rp : double, 3xN
%   Cartesian coordinates of measurement locations in planet--solar--magnetospheric frame in units
%   of planetary radii.
% Psi : double, 1xN
%   Dipole angle relative to perpendicular to pointing into the solar wind in radians, i.e.
%   ``Psi = 0`` is perpendicular to the solar wind direction, ``Psi = pi/2`` is directed into the
%   solar wind, and ``Psi = -pi/2`` is directed along the solar wind velocity vector.
% a, b : double, 3x3
%   Model coefficients fit from Arridge and Eggington (2021).
% c, d : double, 4x4
%   Model coefficients fit from Arridge and Eggington (2021).
% pqrs : double, 4x4
%   Model coefficients fit from Arridge and Eggington (2021). p and r are 3x1, q and s are 4x1.
%
% Returns
% -------
% Bx, By, Bz : double, 1xN
%   Magnetic field contribution from magnetopause surface currents in planet--solar--magnetospheric
%   coordinates.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [Bperpx, Bperpy, Bperpz, Bparax, Bparay, Bparaz] = deal(zeros(1, length(Psi)));
    p = pqrs(1:3, 1);
    q = pqrs(:, 2);
    r = pqrs(1:3, 3);
    s = pqrs(:, 4);

    for i=1:3
        for k=1:3
            prScale = sqrt(1/p(i)^2 + 1/r(k)^2);
            abRot = a(i,k)*cos(Psi) + b(i,k)*cos(2*Psi);
            dBperpx = -prScale * abRot .* exp(xyzPSM_Rp(1,:)*prScale) ...
                .* cos(xyzPSM_Rp(2,:)/p(i)) .* sin(xyzPSM_Rp(3,:)/r(i));
            dBperpy =   1/p(i) * abRot .* exp(xyzPSM_Rp(1,:)*prScale) ...
                .* sin(xyzPSM_Rp(2,:)/p(i)) .* sin(xyzPSM_Rp(3,:)/r(i));
            dBperpz =  -1/r(k) * abRot .* exp(xyzPSM_Rp(1,:)*prScale) ...
                .* cos(xyzPSM_Rp(2,:)/p(i)) .* cos(xyzPSM_Rp(3,:)/r(i));

            Bperpx = Bperpx + dBperpx;
            Bperpy = Bperpy + dBperpy;
            Bperpz = Bperpz + dBperpz;
        end
    end

    for j=1:4
        for l=1:4
            qsScale = sqrt(1/q(j)^2 + 1/s(l)^2);
            cdRot = c(j,l)*sin(Psi) + d(j,l)*sin(2*Psi);
            dBparax = -qsScale * cdRot .* exp(xyzPSM_Rp(1,:)*qsScale) ...
                .* cos(xyzPSM_Rp(2,:)/q(j)) .* cos(xyzPSM_Rp(3,:)/s(l));
            dBparay =   1/q(j) * cdRot .* exp(xyzPSM_Rp(1,:)*qsScale) ...
                .* sin(xyzPSM_Rp(2,:)/q(j)) .* cos(xyzPSM_Rp(3,:)/s(l));
            dBparaz =   1/s(l) * cdRot .* exp(xyzPSM_Rp(1,:)*qsScale) ...
                .* cos(xyzPSM_Rp(2,:)/q(j)) .* sin(xyzPSM_Rp(3,:)/s(l));

            Bparax = Bparax + dBparax;
            Bparay = Bparay + dBparay;
            Bparaz = Bparaz + dBparaz;
        end
    end

    Bx = Bperpx + Bparax;
    By = Bperpy + Bparay;
    Bz = Bperpz + Bparaz;
end
