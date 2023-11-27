function ApplyPlotDefaults(fig, interpreter, font)
% Applies default settings to plot labels.
%
% Parameters
% ----------
% fig : figure
%   Figure object for the plot to adjust.
% interpreter : char, 1xC
%   The default interpreter to use for rendering text labels on plots.
% font : char, 1xD
%   The default font to use for text labels on plots.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    set(fig, 'defaultTextInterpreter', interpreter);  
    set(fig, 'defaultAxesTickLabelInterpreter', interpreter);  
    set(fig, 'defaultLegendInterpreter', interpreter);
    set(fig, 'defaultTextFontName', font);
    set(fig, 'defaultAxesFontName', font);
    set(fig, 'defaultLegendFontName', font);

end
