function [xx, yy, windowName, titleInfo, xInfo, yInfo, fName, xlims, ylims] = HodogramSingle( ...
    moonName, parentName, magModelDescrip, Bvec, SPHOUT, COMPARE_SEUFERT, COMPARE_PHIO)
% Organizes data and labels for plotting a hodogram, showing the projection of the tip of the
% magnetic field vector in a relevant plane, usually the equatorial plane.
%
% Must be plotted with a separate function.
%
% Parameters
% ----------
% moonName : char, 1xC
%   Name of moon for which to evaluate excitation moments.
% parentName : char, 1xD
%   Name of parent planet for the moon for which to evaluate excitation moments.
% magModelDescrip : char, 1xE
%   Text description of the magnetic field model that was evalauted for the input time series.
% Bvec : double, 3xN
%   Magnetic field vectors at each measurement time.
% SPHOUT : bool, default=0
%   Whether to return vectors aligned to spherical coordinate axes (true) or cartesian (false).
% COMPARE_SEUFERT : bool, default=1
%   Whether to plot coordinate directions aligned to the axes presented in Seufert et al. (2011)
%   https://doi.org/10.1016/j.icarus.2011.03.017. Only has an effect when ``SPHOUT`` is true,
%   because Seufert et al. used spherical coordinates for evaluating the excitation spectra of
%   Jupiter's moons.
% COMPARE_PHIO : bool, default=1
%   Whether to plot components in the same orientation as PhiO axes, as in past work, and label the
%   axes with the comparisons (true), or just plot IAU frame axes without comparisons (false).
%
% Returns
% -------
% xx : double, 1xN
%   x axis data.
% yy : double, 1xN
%   y axis data.
% windowName : char, 1xC
%   Name to use for the interactive figure window, which is shown when ``LIVE_PLOTS`` is true.
% titleInfo : char, 1xD
%   Description of plot contents to print at the top of the figure.
% xInfo : char, 1xE
%   Label for x axis of plot.
% yInfo : char, 1xF
%   Label for y axis of plot.
% fName : char, 1xG
%   File name to use for output figures, sans extension.
% xlims : double, 1x2
%   Minimum and maximum values to display on x axis, respectively.
% ylims : double, 1x2
%   Minimum and maximum values to display on y axis, respectively.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~exist('SPHOUT', 'var'); SPHOUT = 0; end
    if ~exist('COMPARE_SEUFERT', 'var'); COMPARE_SEUFERT = 1; end
    if ~exist('COMPARE_PHIO', 'var'); COMPARE_PHIO = 1; end
    % The following are defined in SetPlotDefaults. Do NOT reset them anywhere else.
    global nmTxt
    global bnmTxt
    global mathTxt
    global bmathTxt

    windowName = [moonName ' ' magModelDescrip ' hodogram'];
    if SPHOUT
        if strcmp(parentName,'Saturn')
            BcompMax = max(abs(Bvec(1:2,:)), [], 'all');
            titleInfo = [bnmTxt moonName ' spin-parent plane hodogram, ' magModelDescrip];
            if COMPARE_PHIO
                xx = -Bvec(1,:);
                yy = -Bvec(2,:);
                xInfo = [mathTxt '-B_r ' nmTxt parentName ' SIII (' mathTxt '\approx B_y ' ...
                    nmTxt moonName(1) mathTxt '\phi\Omega' nmTxt ', nT)'];
                yInfo = [mathTxt '-B_\theta ' nmTxt parentName ' SIII (' mathTxt '\approx B_z ' ...
                    nmTxt moonName(1) mathTxt '\phi\Omega' nmTxt ', nT)'];
            else
                xx = Bvec(1,:);
                yy = Bvec(2,:);
                xInfo = [mathTxt 'B_r ' nmTxt parentName ' SIII (nT)'];
                yInfo = [mathTxt 'B_\theta ' nmTxt parentName ' SIII (nT)'];
            end
        else
            BcompMax = max(abs([Bvec(1,:), Bvec(3,:)]));
            titleInfo = [moonName ' equatorial plane hodogram, ' magModelDescrip];
            if COMPARE_SEUFERT
                xx = -Bvec(3,:);
                yy =  Bvec(1,:);
                if COMPARE_PHIO
                    xInfo = [mathTxt '-B_\phi ' nmTxt parentName ' SIII (' mathTxt ...
                        '\approx -B_x ' nmTxt moonName(1) mathTxt '\phi\Omega' nmTxt ', nT)'];
                    yInfo = [mathTxt 'B_r ' nmTxt parentName ' SIII (' mathTxt ...
                        '\approx -B_y ' nmTxt moonName(1) mathTxt '\phi\Omega' nmTxt ', nT)'];
                else
                    xInfo = [mathTxt '-B_\phi ' nmTxt parentName ' SIII (nT)'];
                    yInfo = [mathTxt 'B_r ' nmTxt parentName ' SIII (nT)'];
                end
            else
                if COMPARE_PHIO
                    xx =  Bvec(3,:);
                    yy = -Bvec(1,:);
                    xInfo = [mathTxt 'B_\phi ' nmTxt parentName ' SIII (' mathTxt ...
                        '\approx B_x ' nmTxt moonName(1) mathTxt '\phi\Omega' nmTxt ', nT)'];
                    yInfo = [mathTxt '-B_r ' nmTxt parentName ' SIII (' mathTxt ...
                        '\approx B_y ' nmTxt moonName(1) mathTxt '\phi\Omega' nmTxt ', nT)'];
                else
                    xx = Bvec(1,:);
                    yy = Bvec(3,:);
                    xInfo = [mathTxt 'B_r ' nmTxt parentName ' SIII (nT)'];
                    yInfo = [mathTxt 'B_\phi ' nmTxt parentName ' SIII (nT)'];
                end
            end
        end
    else
        if strcmp(parentName,'Saturn')
            BcompMax = max(abs([Bvec(1,:), Bvec(3,:)]));
            titleInfo = [moonName ' spin-parent plane hodogram, ' magModelDescrip];
            xx = Bvec(1,:);
            yy = Bvec(3,:);
            if COMPARE_PHIO
                xInfo = [mathTxt 'B_x ' nmTxt 'IAU (' mathTxt '\approx B_y ' nmTxt moonName(1) ...
                    mathTxt '\phi\Omega' nmTxt ', nT)'];
                yInfo = [mathTxt 'B_z ' nmTxt 'IAU (' mathTxt '\approx B_z ' nmTxt moonName(1) ...
                    mathTxt '\phi\Omega' nmTxt ', nT)'];
            else
                xInfo = [mathTxt 'B_x ' nmTxt 'IAU (nT)'];
                yInfo = [mathTxt 'B_z ' nmTxt 'IAU (nT)'];
            end
        else
            BcompMax = max(abs(Bvec(1:2,:)), [], 'all');
            titleInfo = [moonName ' equatorial plane hodogram, ' magModelDescrip];
            if COMPARE_PHIO
                xx = Bvec(2,:);
                yy = Bvec(1,:);
                xInfo = [mathTxt 'B_y ' nmTxt 'IAU (' mathTxt '\approx -B_x ' nmTxt moonName(1) ...
                    mathTxt '\phi\Omega' nmTxt ', nT)'];
                yInfo = [mathTxt 'B_x ' nmTxt 'IAU (' mathTxt '\approx B_y ' nmTxt moonName(1) ...
                    mathTxt '\phi\Omega' nmTxt ', nT)'];
            else
                xx = Bvec(1,:);
                yy = Bvec(2,:);
                xInfo = [mathTxt 'B_x ' nmTxt 'IAU (nT)'];
                yInfo = [mathTxt 'B_y ' nmTxt 'IAU (nT)'];
            end
        end
    end
    ylims = [-BcompMax, BcompMax] * 1.05;
    xlims = ylims;

    fName = [moonName 'Hodogram' magModelDescrip];
end