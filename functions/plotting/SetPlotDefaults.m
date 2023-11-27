function [interpreter, font] = SetPlotDefaults
% Sets global Matlab plotting defaults, including some global font switches.
%
% Returns
% -------
% interpreter : char, 1xC
%   The default interpreter to use for rendering text labels on plots.
% font : char, 1xD
%   The default font to use for text labels on plots.
% nmTxt : char, 1xE (global)
%   Font switch to use in LaTeX-formatted labels with "normal" text.
% bnmTxt : char, 1xF (global)
%   Font switch to use in LaTeX-formatted title labels (with bold text).
% mathTxt : char, 1xG (global)
%   Font switch to use in LaTeX-formatted labels with math text.
% bmathTxt : char, 1xH (global)
%   Font switch to use in LaTeX-formatted labels with bolded math text.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    global nmTxt
    global bnmTxt
    global mathTxt
    global bmathTxt
    nmTxt   = '\rm{}';
    bnmTxt  = '\rm\bf';
    mathTxt = '\it{}';
    bmathTxt  = '\bf{}';

    interpreter = 'tex';
    font = 'Times';

end
