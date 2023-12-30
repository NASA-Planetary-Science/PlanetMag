function BD = ICAdecomposition(moonName, parentName, ets, Bvec, magModelDescrip, SPHOUT, ...
    PLOT_DIAGNOSTIC, COMPARE_SEUFERT, COMPARE_PHIO, LIVE_PLOTS, figDir, figXtn)
% Decomposes the input magnetic field vector time series into complex excitation moments using
% :dfn:`Independent Component Analysis (ICA)`.
%
% The methodology applied to invert the imput magnetic field vector time series ``Bvec`` for its
% complex excitation moments (the "independent components" of ICA) is described in detail in
% Hyvarinen and Oja (2000) https://doi.org/10.1016/S0893-6080(00)00026-5. In this method, the
% expected signal components :math:`\mathbf{s} = (s_1, \dots, s_p)` are found from the input
% measurements :math:`\mathbf{x} = (x_1, \dots, x_n)` by inverting the estimated mixing matrix
% :math:`\mathbf{A}` and multiplying :math:`\mathbf{A}^{-1}` by the input time series, where
% :math:`\mathbf{x} = \mathbf{As}`. The signal components :math:`\mathbf{s}` in our case are a list
% of cos and sin waves, representing the real and imaginary parts of :math:`e^{-i\omega_k t}`
% respectively, oscillating with the expected frequencies :math:`\mathbf{f} = (f_1, \dots, f_p)`,
% evaluated at each ephemeris time in the input time series ``ets``. Weighting the signal
% components by the excitation moments :math:`\mathbf{B}_j^e` reproduces the input signal
% :math:`\mathbf{B}_j` when the excitation moments are perfectly resolved. Each signal component
% comprises a column of the signal matrix :math:`\mathbf{X}`. The eigenvectors of
% :math:`\mathbf{X}` are the columns of the weight matrix 
% :math:`\mathbf{W} = (\mathbf{X}^T\mathbf{X})^{-1}`. We find the excitation moments from
%
% .. math::
%   \mathbf{B}_j^e = \mathbf{B}_j\mathbf{XW},
%
% with each element :math:`B_{jk}^e` corresponding to a real (cos) or imaginary (sin) wave with
% frequency :math:`f_k`. The complex excitation moments are constructed by combining the real and
% imaginary parts for each :math:`f_k`, and the reproduced time series is found by taking the real
% part of
%
% .. math::
%   \tilde{B}_j(t) = \sum_k \tilde{B}_{jk}^e e^{-i\omega_k t},
%
% or for the input time series
%
% .. math::
%   \tilde{\mathbf{B}}_j = \mathbf{X}(\mathbf{B}_j^e)^T,
%
% where :math:`\tilde{B}` denotes a complex magnetic field quantity.
%
% Parameters
% ----------
% moonName : char, 1xC
%   Name of moon for which to evaluate excitation moments.
% parentName : char, 1xD
%   Name of parent planet for the moon for which to evaluate excitation moments.
% ets : double, 1xN
%   Ephemeris times associated with input time series magnetic field vectors in TDB seconds
%   relative to J2000.
% Bvec : double, 3xN
%   Magnetic field vectors at each measurement time.
% magModelDescrip : char, 1xE
%   Text description of the magnetic field model that was evalauted for the input time series.
% SPHOUT : bool, default=0
%   Whether to return vectors aligned to spherical coordinate axes (true) or cartesian (false).
% PLOT_DIAGNOSTIC : bool, default=1
%   Whether to plot comparisons between the input and reproduced time series.
% COMPARE_SEUFERT : bool, default=1
%   Whether to plot coordinate directions aligned to the axes presented in Seufert et al. (2011)
%   https://doi.org/10.1016/j.icarus.2011.03.017. Only has an effect when ``SPHOUT`` is true,
%   because Seufert et al. used spherical coordinates for evaluating the excitation spectra of
%   Jupiter's moons.
% COMPARE_PHIO : bool, default=1
%   Whether to plot components in the same orientation as PhiO axes, as in past work, and label the
%   axes with the comparisons (true), or just plot IAU frame axes without comparisons (false).
% LIVE_PLOTS : bool, default=0
%   Whether to load interactive figure windows for plots (true) or print them to disk (false).
% figDir : char, 1xF, default='figures'
%   Whether to load interactive figure windows for plots (true) or print them to disk (false).
% figXtn : char, 1xG, default='pdf'
%   Extension to use for figures, which determines the file type.
%
% Returns
% -------
% BD : struct
%   Contains fields:
%
%       - **f** (`double, 1xP`) -- Frequencies of estimated excitation moments in Hz.
%       - **BexcVec1i** (`double, 1xP`) -- Real part (i for "in-phase") of excitation moments for
%         vector component 1 (x or r) in nT.
%       - **BexcVec1q** (`double, 1xP`) -- Imaginary part (q for "quadrature") of excitation
%         moments for vector component 1 in nT.
%       - **BexcVec1o** (`double`) -- Vector component 1 of uniform background field in nT.
%       - **BexcVec2i** (`double, 1xP`) -- Real part (i for "in-phase") of excitation moments for
%         vector component 2 (y or theta) in nT.
%       - **BexcVec2q** (`double, 1xP`) -- Imaginary part (q for "quadrature") of excitation
%         moments for vector component 2 in nT.
%       - **BexcVec2o** (`double`) -- Vector component 2 of uniform background field in nT.
%       - **BexcVec3i** (`double, 1xP`) -- Real part (i for "in-phase") of excitation moments for
%         vector component 3 (z or phi) in nT.
%       - **BexcVec3q** (`double, 1xP`) -- Imaginary part (q for "quadrature") of excitation
%         moments for vector component 3 in nT.
%       - **BexcVec3o** (`double`) -- Vector component 3 of uniform background field in nT.
%       - **rmse** (`double`) -- Normalized root-mean-squared error calculated from comparing input
%         and reproduced time series data.
%       - **Rsquared** (`double`) -- :math:`R^2` value for the reproduced time series.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~exist('SPHOUT', 'var'); SPHOUT = 0; end
    if ~exist('PLOT_DIAGNOSTIC', 'var'); PLOT_DIAGNOSTIC = 1; end
    if ~exist('COMPARE_SEUFERT', 'var'); COMPARE_SEUFERT = 1; end
    if ~exist('COMPARE_PHIO', 'var'); COMPARE_PHIO = 1; end
    if ~exist('LIVE_PLOTS', 'var'); LIVE_PLOTS = 0; end
    if ~exist('figDir', 'var'); figDir = 'figures'; end
    if ~exist('figXtn', 'var'); figXtn = 'pdf'; end
    % The following are defined in SetPlotDefaults. Do NOT reset them anywhere else.
    global nmTxt
    global bnmTxt
    global mathTxt
    global bmathTxt

    ets=ets(:); % Ensures ets is a row vector
    npts = length(ets);
    ets_ca = ets(round(length(ets)/2));
    % ets_ca = 1228738876.0; % For Triton flyby at AOL 50 degrees
    etMid_day = ets_ca / 86400;

    f = GetExcitations(moonName, etMid_day);

    nFreqs = length(f);
    % Initialize frequency sampling matrix
    X = zeros(npts, 2*nFreqs+1);
    % Fill matrix with sampling for real and imaginary parts for each freq
    for i=1:nFreqs
        X(:, (2*i-1):(2*i)) = [cos(2*pi*f(i)*ets), sin(2*pi*f(i)*ets)];
    end
    % Also fill last column with 1s to get static background field
    X(:, end) = ones(size(ets));

    % Calculate weights, i.e. eigenvectors of X as columns
    W = (X'*X)^-1;
    % Multiply the time series, weighted by each signal component, by the weights to get the
    % amplitudes of excitation
    BexcVec1Est = Bvec(1,:)*X*W;
    BexcVec2Est = Bvec(2,:)*X*W;
    BexcVec3Est = Bvec(3,:)*X*W;
    Bvec1est = X*BexcVec1Est';
    Bvec2est = X*BexcVec2Est';
    Bvec3est = X*BexcVec3Est';

    % Split up the excitation amplitudes into real (i for "in-phase") and imaginary (q for
    % "quadrature") parts
    BexcVec1i = BexcVec1Est(1:2:end-2);
    BexcVec1q = BexcVec1Est(2:2:end-1);
    BexcVec1o = BexcVec1Est(end);
    BexcVec2i = BexcVec2Est(1:2:end-2);
    BexcVec2q = BexcVec2Est(2:2:end-1);
    BexcVec2o = BexcVec2Est(end);
    BexcVec3i = BexcVec3Est(1:2:end-2);
    BexcVec3q = BexcVec3Est(2:2:end-1);
    BexcVec3o = BexcVec3Est(end);

    % Calculate error in reproduction
    Bvec1diff = Bvec(1,:) - Bvec1est';
    Bvec2diff = Bvec(2,:) - Bvec2est';
    Bvec3diff = Bvec(3,:) - Bvec3est';
    BD.rmse = sqrt(sum(Bvec1diff.^2 + Bvec2diff.^2 + Bvec3diff.^2)) / 3/npts;
    dataMean = sum(Bvec, 'all') / 3/npts;
    rmseData = sqrt(sum((Bvec(1,:) - dataMean).^2 + (Bvec(2,:) - dataMean).^2 + (Bvec(3,:) ...
        - dataMean).^2)) / 3/npts;
    BD.Rsquared = 1 - (BD.rmse/rmseData)^2;

    % Fill output struct
    BD.BexcVec1i = BexcVec1i;
    BD.BexcVec1q = BexcVec1q;
    BD.BexcVec1o = BexcVec1o;
    BD.BexcVec2i = BexcVec2i;
    BD.BexcVec2q = BexcVec2q;
    BD.BexcVec2o = BexcVec2o;
    BD.BexcVec3i = BexcVec3i;
    BD.BexcVec3q = BexcVec3q;
    BD.BexcVec3o = BexcVec3o;

    BD.f = f;

    BexcVec1abs = sqrt(BexcVec1i.^2 + BexcVec1q.^2);
    BexcVec2abs = sqrt(BexcVec2i.^2 + BexcVec2q.^2);
    BexcVec3abs = sqrt(BexcVec3i.^2 + BexcVec3q.^2);
    BexcAbs2 = BexcVec1abs.^2 + BexcVec2abs.^2 + BexcVec3abs.^2;
    [~, iMax] = max(BexcAbs2);
    fMax = f(iMax);
    [~, iMaxComp] = max([BexcVec1abs(iMax), BexcVec2abs(iMax), BexcVec3abs(iMax)]);

    disp(' ')
    disp(['R^2 value: ' num2str(BD.Rsquared)])
    disp(' ')
    disp('Driving Field Waves...')
    for i = 1:nFreqs
        disp([num2str(i), '. T = ', num2str(1/f(i)/3600,'%0.2f'), ' hr'])
        % disp(['  2*pi*f*t = ', num2str(2*pi*f(i)*ets_ca*180/pi,'%0.2f'), ' deg'])
        if SPHOUT
            Bv1name = 'Br'; Bv2name = 'Bth'; Bv3name = 'Bphi';
            Pv1name = 'Pr'; Pv2name = 'Pth'; Pv3name = 'Pphi';
            Bv1lbl = 'B_r'; Bv2lbl = 'B_\theta'; Bv3lbl = 'B_\phi';
        else
            Bv1name = 'Bx'; Bv2name = 'By'; Bv3name = 'Bz';
            Pv1name = 'Px'; Pv2name = 'Py'; Pv3name = 'Pz';
            Bv1lbl = 'B_x'; Bv2lbl = 'B_y'; Bv3lbl = 'B_z';
        end
        disp(['  ' Bv1name ' = ' num2str(BexcVec1abs(i),'%0.2f') ' nT, ' Bv2name ' = ' ...
            num2str(BexcVec2abs(i),'%0.2f') ' nT, ' Bv3name ' = ' ...
            num2str(BexcVec3abs(i),'%0.2f') ' nT'])
        disp(['  ' Pv1name ' = ' num2str(atan2(BexcVec1q(i),BexcVec1i(i))*180/pi,'%0.2f') ...
            ' deg, ' Pv2name ' = ' num2str(atan2(BexcVec2q(i),BexcVec2i(i))*180/pi,'%0.2f') ...
            ' deg, ' Pv3name ' = ' num2str(atan2(BexcVec3q(i),BexcVec3i(i))*180/pi,'%0.2f') ...
            ' deg'])
    end
    disp(' ')

    assignin('base','BexcVec1Est',BexcVec1Est)
    assignin('base','BexcVec2Est',BexcVec2Est)
    assignin('base','BexcVec3Est',BexcVec3Est)

    if PLOT_DIAGNOSTIC

        etsRel_h = (ets-ets(1))/3600;

        cosSyn = X(:, 2*(find(f==fMax))-1);
        sinSyn = X(:, 2*(find(f==fMax)));
        Bvec1max = BexcVec1o + BexcVec1i(f==fMax)*cosSyn + BexcVec1q(f==fMax)*sinSyn;
        Bvec2max = BexcVec2o + BexcVec2i(f==fMax)*cosSyn + BexcVec2q(f==fMax)*sinSyn;
        Bvec3max = BexcVec3o + BexcVec3i(f==fMax)*cosSyn + BexcVec3q(f==fMax)*sinSyn;

        switch(iMaxComp)
            case 1
                BvecCompMax = Bvec1max;
            case 2
                BvecCompMax = Bvec2max;
            case 3
                BvecCompMax = Bvec3max;
        end

        figNumBase = 10;
        windowName = [magModelDescrip ' model vs. reconstruction, first 200h'];
        titleInfo = [bnmTxt moonName ' ' magModelDescrip ' model vs. reconstruction, first 200h'];
        xx = etsRel_h;
        yy = [Bvec(1,:); Bvec(2,:); Bvec(3,:); Bvec1est'; Bvec2est'; Bvec3est'; Bvec1max'; ...
            Bvec2max'; Bvec3max'];
        cfmt = {'b', 'r', 'g', '--k', '--k', '--k', 'c', 'm', 'y'};
        legendStrings = [string(Bv1lbl), string(Bv2lbl), string(Bv3lbl), "Reproduced", ...
            "Largest only"];
        xInfo = 'Time (hr)';
        yInfo = 'Magnetic Field (nT)';
        xlims = [0 200];
        fName = [moonName 'BVecReconVs' magModelDescrip];
        PlotGeneric(xx, yy, legendStrings, windowName, titleInfo, xInfo, yInfo, fName, ...
            figDir, figXtn, LIVE_PLOTS, figNumBase + 1, 'linear', 'linear', xlims, 'auto', ...
            cfmt, 1, 1);

        windowName = [magModelDescrip ' model vs. reconstruction diff, first 200h'];
        titleInfo = [bnmTxt moonName ' ' magModelDescrip ' model vs. reconstruction diff, ' ...
            'first 200h'];
        yy = [Bvec1diff; Bvec2diff; Bvec3diff; Bvec(iMaxComp,:) - BvecCompMax'];
        cfmt = {'b', 'r', 'g', 'm'};
        legendStrings = [string(Bv1lbl), string(Bv2lbl), string(Bv3lbl), "Largest only"];
        yInfo = 'Magnetic Field Error (nT)';
        fName = [moonName 'DeltaBreconVs' magModelDescrip];
        PlotGeneric(xx, yy, legendStrings, windowName, titleInfo, xInfo, yInfo, fName, ...
            figDir, figXtn, LIVE_PLOTS, figNumBase + 2, 'linear', 'linear', xlims, 'auto', ...
            cfmt, 1, 1);

        %% Hodograms
        [xx, yy, windowName, titleInfo, xInfo, yInfo, fName, xlims, ylims] = HodogramSingle( ...
            moonName, parentName, magModelDescrip, Bvec, SPHOUT, COMPARE_SEUFERT, COMPARE_PHIO);
        PlotGeneric(xx, yy, [], windowName, titleInfo, xInfo, yInfo, fName, figDir, figXtn, ...
            LIVE_PLOTS, figNumBase + 3, 'linear', 'linear', xlims, ylims, {'b'}, 1, 1);
    end

end
