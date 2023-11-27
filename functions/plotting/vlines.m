function vlines(x, color)
% Plot vertical lines at specific x locations with the named color.
%
% Vertical lines are plotted on the currently focused axes object.
%
% Parameters
% ----------
% x : double, 1xN
%   x coordinates at which to plot vertical lines in data units.
% color : char, 1x7
%   Object interpretable by Matlab as a color (e.g. HTML color string ``'#FF0000'``).

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    yLim = get(gca,'ylim');
    for i=1:length(x)
        plot([x(i), x(i)], yLim, 'Color', color);
    end
end
