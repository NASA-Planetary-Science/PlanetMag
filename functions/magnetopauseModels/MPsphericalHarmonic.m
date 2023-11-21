function [Br, Bth, Bphi] = MPsphericalHarmonic(r_Rss, thMP, phiMP, Gnm, Nmax)
% Evaluate external-source spherical harmonics for magnetopause surface currents, i.e. all hnm = 0.
%
% Parameters
% ----------
% r_Rss : double, 1xN
%   Radial distance from planet center of mass in units of the sub-solar magnetopause standoff
%   distance.
% thMP : double, 1xN
%   Colatitude in solar--magnetospheric coordinates. The :math:`\hat{z}` axis of this coordinate
%   system is derived from :math:`\hat{x}\times\hat{y}`, where :math:`\hat{x}` the along the
%   planet--Sun direction and :math:`\hat{y}` is along :math:`\mathbf{m}\times\hat{x}`, with
%   :math:`\mathbf{m}` the instantaneous magnetic dipole moment vector for the body.
% phiMP : double, 1xN
%   Azimuthal angle in solar-magnetospheric coordinates, i.e. the angle between the planet--Sun
%   direction and the position vector's projection into the dipole equatorial plane.
% Gnm : double, (Nmax)x(Nmax+1)
%   Coefficients for external spherical harmonics in nT at body center.
% Nmax : int
%   Maximum degree to which to limit spherical harmonic models. Only has an effect if the value
%   passed is less than the lesser of the model maximum degree.
%
% Returns
% -------
% Br, Bth, Bphi : double, 1xN
%   External-source magnetic field vector components ailgned to spherical coordinates.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [Br, Bth, Bphi] = deal(zeros(1,length(r_Rss)));
    for n=1:min([10, Nmax])

        % Handle m=0 separately to avoid unneccessary calculations for Bphi, which is always 0 for
        % m=0
        Pn0 =   LegendreS(n, 0, thMP);
        dPn0 = dLegendreS(n, 0, thMP);

        dBr = n*(r_Rss).^(n-1) .* squeeze(Gnm(n,1,:))' .*  Pn0;
        dBth =  (r_Rss).^(n-1) .* squeeze(Gnm(n,1,:))' .* dPn0;
        Br = Br + dBr;
        Bth = Bth + dBth;

        for m=1:n

            % Get Schmidt semi-normalized Legendre functions
            thMPadj = thMP;
            thMPadj(thMPadj == 0) = eps;
            Pnm =   LegendreS(n, m, thMPadj);
            dPnm = dLegendreS(n, m, thMP);

            dBr =    n*(r_Rss).^(n-1) .* squeeze(Gnm(n,m+1,:))' .*  Pnm .* cos(m*phiMP);
            dBth =     (r_Rss).^(n-1) .* squeeze(Gnm(n,m+1,:))' .* dPnm .* cos(m*phiMP);
            dBphi = -m*(r_Rss).^(n-1) .* squeeze(Gnm(n,m+1,:))' .*  Pnm ./ sin(thMPadj) ...
                .* sin(m*phiMP);

            Br = Br + dBr;
            Bth = Bth + dBth;
            Bphi = Bphi + dBphi;

        end
    end
end
