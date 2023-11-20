function BD = ICAdecomposition(ets, moonName, parentName, Bvec, descrip, SPHOUT, ...
    PLOT_DIAGNOSTIC, COMPARE_SEUFERT, LIVE_PLOTS)
% Decomposes the input magnetic field vector time series into complex excitation moments using
% :dfn:`Independent Component Analysis (ICA)`.
%
% The methodology applied to invert the imput magnetic field vector time series ``Bvec`` for its
% complex excitation moments (the independent components of ICA) is described in detail in Hyvarinen
% and Oja (2000) https://doi.org/10.1016/S0893-6080(00)00026-5. In this method, the expected signal
% components :math:`\mathbf{s} = (s_1, \dots, s_p)`, are found from the input measurements :math:`\mathbf{x} = (x_1, \dots, x_n)` by inverting the estimated mixing matrix :math:`\mathbf{A}` and multiplying :math:`\mathbf{A}^{-1}` by the input time series, where :math:`\mathbf{x} = \mathbf{As}`. The signal components :math:`\mathbf{s}_j` for vector component :math:`j` in our case are a list of cos and sin waves, representing the real and imaginary parts of :math:`e^{-i\omega t}`, respectively, oscillating with the expected frequencies :math:`\mathbf{f} = (f_1, \dots, f_p)`, evaluated at each ephemeris time in the input time series ``ets``, and weighted by the excitation moments :math:`\mathbf{B}_j^e`. Each unweighted signal component comprises a column of :math:`\mathbf{X}`. The eigenvectors of :math:`\mathbf{X}` are the columns of the weight matrix :math:`\mathbf{W} = (\mathbf{X}^T \mathbf{X})^{-1}`. We find the excitation moments from :math:`\mathbf{B}_j^e = \mathbf{X}(\mathbf{B}_j\mathbf{XW})^T`, because each row of this product describes the
% oscillation frequencies :math:`f_k` are each used to create "principal
% components", the sinusoidal oscillations associated
% with each frequency (both cos and sin, for the ). The principal components are used to construct a
% matrix, X, where each principal component is a row in X. The inverse of
% X^T * X is then calculated, which yields the eigenvectors as columns.
% Multiplying (Bi*X) by the column-eigenvector matrix then weights the 
% principal components according to their degree of representation in the
% input series (where Bi is one of Bx/By/Bz or Br/Btheta/Bphi etc). This calculation is made
% for each B component to get the excitation moments for each.
%
% Parameters
% ----------
% SPHOUT : bool, default=0
%   Whether to return vectors aligned to spherical coordinate axes (true) or cartesian (false).
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

    if ~exist('SPHOUT', 'var'); SPHOUT = 0; end
    if ~exist('PLOT_DIAGNOSTIC', 'var'); PLOT_DIAGNOSTIC = 1; end
    if ~exist('COMPARE_SEUFERT', 'var'); COMPARE_SEUFERT = 0; end
    if ~exist('LIVE_PLOTS', 'var'); LIVE_PLOTS = 0; end

    ets=ets(:); % Ensures ets is a row vector
    npts = length(ets);
    ets_ca = ets(round(length(ets)/2));
    % ets_ca = 1228738876.0; % For Triton flyby at AOL 50 degrees
    etMid_day = ets_ca / 86400;

    [f, fSyn] = GetExcitations(moonName, etMid_day);

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
    Bdvec1Est = Bvec(1,:)*X*W;
    Bdvec2Est = Bvec(2,:)*X*W;
    Bdvec3Est = Bvec(3,:)*X*W;
    Bvec1est = X*Bdvec1Est';
    Bvec2est = X*Bdvec2Est';
    Bvec3est = X*Bdvec3Est';

    % Split up the excitation amplitudes into real (i for "in-phase") and imaginary (q for
    % "quadrature") parts
    Bdvec1i = Bdvec1Est(1:2:end-2);
    Bdvec1q = Bdvec1Est(2:2:end-1);
    Bdvec1o = Bdvec1Est(end);
    Bdvec2i = Bdvec2Est(1:2:end-2);
    Bdvec2q = Bdvec2Est(2:2:end-1);
    Bdvec2o = Bdvec2Est(end);
    Bdvec3i = Bdvec3Est(1:2:end-2);
    Bdvec3q = Bdvec3Est(2:2:end-1);
    Bdvec3o = Bdvec3Est(end);

    % Fill output struct
    BD.Bdvec1i = Bdvec1i;
    BD.Bdvec1q = Bdvec1q;
    BD.Bdvec1o = Bdvec1o;
    BD.Bdvec2i = Bdvec2i;
    BD.Bdvec2q = Bdvec2q;
    BD.Bdvec2o = Bdvec2o;
    BD.Bdvec3i = Bdvec3i;
    BD.Bdvec3q = Bdvec3q;
    BD.Bdvec3o = Bdvec3o;

    BD.f = f;

    disp(' ')
    disp('Driving Field Waves...')
    for i = 1:nFreqs
        disp([num2str(i), '. T = ', num2str(1/f(i)/3600,'%0.2f'), ' hr'])
        Bdvec1 = sqrt(Bdvec1i(i)^2+Bdvec1q(i)^2);
        Bdvec2 = sqrt(Bdvec2i(i)^2+Bdvec2q(i)^2);
        Bdvec3 = sqrt(Bdvec3i(i)^2+Bdvec3q(i)^2);
        % disp(['  2*pi*f*t = ', num2str(2*pi*f(i)*ets_ca*180/pi,'%0.2f'), ' deg'])
        if SPHOUT
            Bv1name = 'Br'; Bv2name = 'Bth'; Bv3name = 'Bphi';
            Pv1name = 'Pr'; Pv2name = 'Pth'; Pv3name = 'Pphi';
        else
            Bv1name = 'Bx'; Bv2name = 'By'; Bv3name = 'Bz';
            Pv1name = 'Px'; Pv2name = 'Py'; Pv3name = 'Pz';
        end
        disp(['  ' Bv1name ' = ' num2str(Bdvec1,'%0.2f') ' nT, ' Bv2name ' = ' ...
            num2str(Bdvec2,'%0.2f') ' nT, ' Bv3name ' = ' num2str(Bdvec3,'%0.2f') ' nT'])
        disp(['  ' Pv1name ' = ' num2str(atan2(Bdvec1q(i),Bdvec1i(i))*180/pi,'%0.2f') ' deg, ' ...
            Pv2name ' = ' num2str(atan2(Bdvec2q(i),Bdvec2i(i))*180/pi,'%0.2f') ' deg, ' Pv3name ...
            ' = ' num2str(atan2(Bdvec3q(i),Bdvec3i(i))*180/pi,'%0.2f') ' deg'])
    end
    disp(' ')

    assignin('base','Bdvec1Est',Bdvec1Est)
    assignin('base','Bdvec2Est',Bdvec2Est)
    assignin('base','Bdvec3Est',Bdvec3Est)

    if PLOT_DIAGNOSTIC

    %     figure
    %     hold on; grid on; box on;
    %     for i=1:MagFreqs
    %         plot([1/f(i)/3600,1/f(i)/3600],[1e-8,abs(Bdxi(i))],'b','linewidth',2)
    %         plot([1/f(i)/3600,1/f(i)/3600],[1e-8,abs(Bdxq(i))],'r','linewidth',2)
    %     end
    %     ylabel('Magnetic Field (nT)'); xlabel('Period (hr)');
    %     set(gca,'fontsize',12); set(gca, 'YScale', 'log'); set(gca, 'XScale', 'log');
    %
    %     figure
    %     hold on; grid on; box on;
    %     for i=1:MagFreqs
    %         plot([1/f(i)/3600,1/f(i)/3600],[1e-8,abs(Bdyi(i))],'b','linewidth',2)
    %         plot([1/f(i)/3600,1/f(i)/3600],[1e-8,abs(Bdyq(i))],'r','linewidth',2)
    %     end
    %     ylabel('Magnetic Field (nT)'); xlabel('Period (hr)');
    %     set(gca,'fontsize',12); set(gca, 'YScale', 'log'); set(gca, 'XScale', 'log');

        etsRel_h = (ets-ets(1))/3600;

        if ~strcmp(parentName,'Saturn')
            cosSyn = X(:, 2*(find(f==fSyn))-1);
            sinSyn = X(:, 2*(find(f==fSyn)));
            Bvec1Syn = Bdvec1o + Bdvec1i(f==fSyn)*cosSyn + Bdvec1q(f==fSyn)*sinSyn;
            Bvec2Syn = Bdvec2o + Bdvec2i(f==fSyn)*cosSyn + Bdvec2q(f==fSyn)*sinSyn;
            Bvec3Syn = Bdvec3o + Bdvec3i(f==fSyn)*cosSyn + Bdvec3q(f==fSyn)*sinSyn;
        end

        figure; hold on; box on; grid on;
        set(gcf,'Name', [descrip ' model vs. reconstruction, first 200h']);
        title([moonName ' ' descrip ' model vs. reconstruction, first 200h'])
        Bvec1in = plot(etsRel_h, Bvec(1,:), 'b');
        Bvec2in = plot(etsRel_h, Bvec(2,:), 'r');
        Bvec3in = plot(etsRel_h, Bvec(3,:), 'g');
        Bvec1out = plot(etsRel_h, Bvec1est, '--k');
        Bvec2out = plot(etsRel_h, Bvec2est, '--k');
        Bvec3out = plot(etsRel_h, Bvec3est, '--k');
        if ~strcmp(parentName,'Saturn')
            Bvec1synOnly = plot(etsRel_h, Bvec1Syn, 'm');
            Bvec2synOnly = plot(etsRel_h, Bvec2Syn, 'm');
            Bvec3synOnly = plot(etsRel_h, Bvec3Syn, 'm');
            legend([Bvec1in,Bvec2in,Bvec3in, Bvec1out, Bvec1synOnly], Bv1name, Bv2name, Bv3name, ...
                'Reproduced', 'Synodic only')
        else
            legend([Bvec1in,Bvec2in,Bvec3in, Bvec3out], Bv1name, Bv2name, Bv3name, 'Reproduced')
        end
        xlabel('Time (hr)')
        ylabel('Magnetic Field (nT)')
        set(gca,'fontsize',16)
        xlim([0 200])

        figure; hold on; box on; grid on;
        set(gcf,'Name', [descrip ' model vs. reconstruction diff, first 200h']);
        title([moonName ' ' descrip ' model vs. reconstruction diff, first 200h'])
        Bvec1diff = plot(etsRel_h, Bvec(1,:)-Bvec1est', 'b');
        Bvec2diff = plot(etsRel_h, Bvec(2,:)-Bvec2est', 'r');
        Bvec3diff = plot(etsRel_h, Bvec(3,:)-Bvec3est', 'g');
        if ~strcmp(parentName,'Saturn')
            Bvec1syndiff = plot(etsRel_h, Bvec(1,:)-Bvec1Syn', 'm');
            legend([Bvec1diff,Bvec2diff,Bvec3diff, Bvec1syndiff], Bv1name, Bv2name, Bv3name,
                [Bv1name 'synodic only'])
        else
            legend([Bvec1diff,Bvec2diff,Bvec3diff], Bv1name, Bv2name, Bv3name)
        end
        xlabel('Time (hr)')
        ylabel('Magnetic Field Error (nT)')
        set(gca,'fontsize',16)
        xlim([0 200])

        figure; box on;
        set(gcf,'Name', [moonName ' ' descrip ' hodogram']);
        if SPHOUT
            if strcmp(parentName,'Saturn')
                plot(-Bvec(1,:), -Bvec(2,:), 'b')
                xlabel(['-B_r ' parentName ' SIII (\approx B_y ' moonName(1) '\phi\Omega, nT)'])
                ylabel(['-B_\theta ' parentName ' SIII (\approx B_z ' moonName(1) ...
                    '\phi\Omega, nT)'])
                title([moonName ' spin-parent plane hodogram, ' descrip])
            else
                if COMPARE_SEUFERT
                    plot(-Bvec(3,:), Bvec(1,:), 'b')
                    xlabel(['-B_\phi ' parentName ' SIII (\approx -B_x ' moonName(1) ...
                        '\phi\Omega, nT)'])
                    ylabel(['B_r ' parentName ' SIII (\approx -B_y ' moonName(1) '\phi\Omega, nT)'])
                else
                    plot(Bvec(3,:), -Bvec(1,:), 'b')
                    xlabel(['B_\phi ' parentName ' SIII (\approx B_x ' moonName(1) ...
                        '\phi\Omega, nT)'])
                    ylabel(['-B_r ' parentName ' SIII (\approx B_y ' moonName(1) '\phi\Omega, nT)'])
                end
                title([moonName ' equatorial plane hodogram, ' descrip])
            end
        else
            if strcmp(parentName,'Saturn')
                plot(Bvec(1,:), Bvec(3,:), 'b')
                xlabel(['B_x IAU (\approx B_y ' moonName(1) '\phi\Omega, nT)'])
                ylabel(['B_z IAU (\approx B_z ' moonName(1) '\phi\Omega, nT)'])
                title([moonName ' spin-parent plane hodogram, ' descrip])
            else
                plot(Bvec(2,:), Bvec(1,:), 'b')
                ylabel(['B_x IAU (\approx B_y ' moonName(1) '\phi\Omega, nT)'])
                xlabel(['B_y IAU (\approx -B_x ' moonName(1) '\phi\Omega, nT)'])
                title([moonName ' equatorial plane hodogram, ' descrip])
            end
        end
        ylim([-max(abs(ylim())), max(abs(ylim()))])
        xlim(ylim())
        grid on;
        set(gca,'fontsize',16)
    end

end
