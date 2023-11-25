function rMP_km = GetMPsurfS1997(Rss_km, xi, thDSZ)
% Determine paraboloid magnetopause surface shape based on Shue et al (1997) model.
%
% Determines shape of magnetopause screening surface current shape from the model detailed by Shue
% et al. (1997) https://doi.org/10.1029/97JA00196. This is a simple paraboloid with an input
% sub-solar magnetopause standoff distance.
%
% Parameters
% ----------
% Rss_km : double, 1xN
%   Sub-solar point magnetopause standoff distance in km.
% xi : double, 1x1 or 1xN
%   Exponent to scale paraboloid shape.
% thDSZ : double, 1xN
%   Colatitude in dipole--solar--zenith coordinates in radians, i.e. the angle between the
%   planet--Sun direction and the evaluation point location vector.
%
% Returns
% -------
% rMP_km : double, 1xN
%   Distance from body center of mass to magnetopause surface for each evaluation point.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    rMP_km = Rss_km .* (2 ./ (1 + cos(thDSZ))).^xi;
end
