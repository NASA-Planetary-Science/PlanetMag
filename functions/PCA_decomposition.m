function BD = PCA_decomposition(ets,body,parentName,Bvec, descrip, SPHOUT, ...
    PLOT_DIAGNOSTIC, COMPARE_SEUFERT)
% Evaluates the excitation moments of the input time series Bx, By, and Bz
% by the method of Principal Component Analysis (PCA).
% In this method, the expected oscillation frequencies are each used to
% create "principal components", the sinusoidal oscillations associated
% with each frequency (both cos and sin, for the real and imaginary part of
% e^-iωt, respectively). The principal components are used to construct a
% matrix, X, where each principal component is a row in X. The inverse of
% X^T * X is then calculated, which yields the eigenvectors as columns.
% Multiplying (Bi*X) by the column-eigenvector matrix then weights the 
% principal components according to their degree of representation in the
% input series (where Bi is one of Bx/By/Bz or Br/Btheta/Bphi etc). This calculation is made
% for each B component to get the excitation moments for each.

if ~exist('PLOT_DIAGNOSTIC', 'var')
    PLOT_DIAGNOSTIC = 1;
end
if ~exist('SPHOUT', 'var')
    SPHOUT = 0;
end
if ~exist('COMPARE_SEUFERT', 'var')
    COMPARE_SEUFERT = 0;
end

ets=ets(:); % Ensures ets is a row vector
npts = length(ets);
ets_ca = ets(round(length(ets)/2));
etMid_day = ets_ca / 86400;
%ets_ca = 1228738876.0; % for Triton flyby at AOL 50  degrees

[f, fSyn] = GetExcitations(body, etMid_day);

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
% Multiply the time series, weighted by each principal component, by the
% weights to get the amplitudes of excitation
Bdvec1Est = Bvec(1,:)*X*W;
Bdvec2Est = Bvec(2,:)*X*W;
Bdvec3Est = Bvec(3,:)*X*W;
Bvec1est = X*Bdvec1Est';
Bvec2est = X*Bdvec2Est';
Bvec3est = X*Bdvec3Est';

% Split up the excitation amplitudes into real (i for "in-phase") and
% imaginary (q for "quadrature") parts
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
    Bdvec1 = sqrt(Bdvec1i(i)^2+Bdvec1q(i)^2); Bdvec2 = sqrt(Bdvec2i(i)^2+Bdvec2q(i)^2); Bdvec3 = sqrt(Bdvec3i(i)^2+Bdvec3q(i)^2);
    %disp(['  2*pi*f*t = ', num2str(2*pi*f(i)*ets_ca*180/pi,'%0.2f'), ' deg'])
    if SPHOUT
        Bv1name = 'Br'; Bv2name = 'Bth'; Bv3name = 'Bphi'; Pv1name = 'Pr'; Pv2name = 'Pth'; Pv3name = 'Pphi';
    else
        Bv1name = 'Bx'; Bv2name = 'By'; Bv3name = 'Bz'; Pv1name = 'Px'; Pv2name = 'Py'; Pv3name = 'Pz';
    end
    disp(['  ' Bv1name ' = ' num2str(Bdvec1,'%0.2f') ' nT, ' Bv2name ' = ' num2str(Bdvec2,'%0.2f') ' nT, ' Bv3name ' = ' num2str(Bdvec3,'%0.2f') ' nT'])
    disp(['  ' Pv1name ' = ' num2str(atan2(Bdvec1q(i),Bdvec1i(i))*180/pi,'%0.2f') ' deg, ' Pv2name ' = ' num2str(atan2(Bdvec2q(i),Bdvec2i(i))*180/pi,'%0.2f') ' deg, ' Pv3name ' = ' num2str(atan2(Bdvec3q(i),Bdvec3i(i))*180/pi,'%0.2f') ' deg'])
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
%     ylabel('Magnetic Field (nT)'); xlabel('Period (hr)'); set(gca,'fontsize',12);set(gca, 'YScale', 'log');set(gca, 'XScale', 'log')
%     
%     figure
%     hold on; grid on; box on;
%     for i=1:MagFreqs
%         plot([1/f(i)/3600,1/f(i)/3600],[1e-8,abs(Bdyi(i))],'b','linewidth',2)
%         plot([1/f(i)/3600,1/f(i)/3600],[1e-8,abs(Bdyq(i))],'r','linewidth',2)
%     end
%     ylabel('Magnetic Field (nT)'); xlabel('Period (hr)'); set(gca,'fontsize',12);set(gca, 'YScale', 'log');set(gca, 'XScale', 'log')

    etsRel_h = (ets-ets(1))/3600;
    bodyname = [body(1) lower(body(2:end))];
    
    if ~strcmp(parentName,'Saturn')
        cosSyn = X(:, 2*(find(f==fSyn))-1);
        sinSyn = X(:, 2*(find(f==fSyn)));
        Bvec1Syn = Bdvec1o + Bdvec1i(f==fSyn)*cosSyn + Bdvec1q(f==fSyn)*sinSyn;
        Bvec2Syn = Bdvec2o + Bdvec2i(f==fSyn)*cosSyn + Bdvec2q(f==fSyn)*sinSyn;
        Bvec3Syn = Bdvec3o + Bdvec3i(f==fSyn)*cosSyn + Bdvec3q(f==fSyn)*sinSyn;
    end

    figure; hold on; box on; grid on;
    set(gcf,'Name', [descrip ' model vs. reconstruction, first 200h']);
    title([bodyname ' ' descrip ' model vs. reconstruction, first 200h'])
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
        legend([Bvec1in,Bvec2in,Bvec3in, Bvec1out, Bvec1synOnly], Bv1name, Bv2name, Bv3name, 'Reproduced', 'Synodic only')
    else
        legend([Bvec1in,Bvec2in,Bvec3in, Bvec3out], Bv1name, Bv2name, Bv3name, 'Reproduced')
    end
    xlabel('Time (hr)')
    ylabel('Magnetic Field (nT)')
    set(gca,'fontsize',16)
    xlim([0 200])
    
    figure; hold on; box on; grid on;
    set(gcf,'Name', [descrip ' model vs. reconstruction diff, first 200h']);
    title([bodyname ' ' descrip ' model vs. reconstruction diff, first 200h'])
    Bvec1diff = plot(etsRel_h, Bvec(1,:)-Bvec1est', 'b');
    Bvec2diff = plot(etsRel_h, Bvec(2,:)-Bvec2est', 'r');
    Bvec3diff = plot(etsRel_h, Bvec(3,:)-Bvec3est', 'g');
    if ~strcmp(parentName,'Saturn')
        Bvec1syndiff = plot(etsRel_h, Bvec(1,:)-Bvec1Syn', 'm');
        legend([Bvec1diff,Bvec2diff,Bvec3diff, Bvec1syndiff], Bv1name, Bv2name, Bv3name, [Bv1name 'synodic only'])
    else
        legend([Bvec1diff,Bvec2diff,Bvec3diff], Bv1name, Bv2name, Bv3name)
    end
    xlabel('Time (hr)')
    ylabel('Magnetic Field Error (nT)')
    set(gca,'fontsize',16)
    xlim([0 200])
    
    figure; box on;
    set(gcf,'Name', [bodyname ' ' descrip ' hodogram']);
    if SPHOUT
        if strcmp(parentName,'Saturn')
            plot(-Bvec(1,:), -Bvec(2,:), 'b')
            xlabel(['-B_r ' parentName ' SIII (\approx B_y ' bodyname(1) '\phi\Omega, nT)'])
            ylabel(['-B_\theta ' parentName ' SIII (\approx B_z ' bodyname(1) '\phi\Omega, nT)'])
            title([bodyname ' spin-parent plane hodogram, ' descrip])
        else
            if COMPARE_SEUFERT
                plot(-Bvec(3,:), Bvec(1,:), 'b')
                xlabel(['-B_\phi ' parentName ' SIII (\approx -B_x ' bodyname(1) '\phi\Omega, nT)'])
                ylabel(['B_r ' parentName ' SIII (\approx -B_y ' bodyname(1) '\phi\Omega, nT)'])
            else
                plot(Bvec(3,:), -Bvec(1,:), 'b')
                xlabel(['B_\phi ' parentName ' SIII (\approx B_x ' bodyname(1) '\phi\Omega, nT)'])
                ylabel(['-B_r ' parentName ' SIII (\approx B_y ' bodyname(1) '\phi\Omega, nT)'])
            end
            title([bodyname ' equatorial plane hodogram, ' descrip])
        end
    else
        if strcmp(parentName,'Saturn')
            plot(Bvec(1,:), Bvec(3,:), 'b')
            xlabel(['B_x IAU (\approx B_y ' bodyname(1) '\phi\Omega, nT)'])
            ylabel(['B_z IAU (\approx B_z ' bodyname(1) '\phi\Omega, nT)'])
            title([bodyname ' spin-parent plane hodogram, ' descrip])
        else
            plot(Bvec(2,:), Bvec(1,:), 'b')
            ylabel(['B_x IAU (\approx B_y ' bodyname(1) '\phi\Omega, nT)'])
            xlabel(['B_y IAU (\approx -B_x ' bodyname(1) '\phi\Omega, nT)'])
            title([bodyname ' equatorial plane hodogram, ' descrip])
        end
    end
    ylim([-max(abs(ylim())), max(abs(ylim()))])
    xlim(ylim())
    grid on;
    set(gca,'fontsize',16)
end

end


%% List of excitation periods to invert over
function [f, fSyn] = GetExcitations(body, etMid_day)
    fJup = 870.536 / 360 / 86400;
    fSat = 810.7939024 / 360 / 86400;
    fUra = abs(-501.1600928) / 360 / 86400;
    fNep = 536.3128492 / 360 / 86400;
    Pnut0Jup_deg = [73.32, 24.62, 283.90, 355.80, 119.90, 229.80, 352.25, 113.35];
    Pnut1Jup_degday = [91472.9, 45137.2, 4850.7, 1191.3, 262.1, 64.3, 2382.6, 6070.0] ...
        / 36525;

    switch body

        %% Jupiter moons 
        case 'IO'
            Pnut0 = Pnut0Jup_deg;
            Pnut1 = Pnut1Jup_degday;
            PnutM = [0, 0, -0.085, -0.022, 0, 0, 0, 0];
            fOrb = 203.4889538 / 360 / 86400; % orbital period
            fOrbBeat = 1 / 3600 / 42.289036514276894252;
            fIoTA = 1 / 3600 / 42.315044531808887029; % true anomaly period
            fOrbMid1 = 1 / 3600 / 42.431381950208163;
            fOrbMid2 = 1 / 3600 / 42.305626942856122;
            fOrbAdj = (203.4889538 + sum(PnutM .* Pnut1 .* cosd(Pnut0 + Pnut1*etMid_day))) / 360 / 86400;
            fSyn = fJup - fOrb;
            fSynAdj = fJup - fOrbAdj;

            f = [fSyn, fOrb, fIoTA, fOrbMid1, fOrbMid2];

            f = [f, 2*fSyn, ...
                    3*fSyn, ...
                    4*fSyn, ...
                    5*fSyn, ...
                    6*fSyn, ...
                    fSyn - fIoTA, ... % 1st harmonic beats
                    fSyn - fOrb, ...
                    fSyn + fOrb, ...
                    fSyn + fIoTA, ...
                    2*fSyn - fIoTA, ... % 2nd harmonic beats
                    2*fSyn + fIoTA, ...
                    3*fSyn - fIoTA, ... % 3rd harmonic beats
                    3*fSyn + fIoTA, ...
                    fSynAdj - fIoTA, ... % 1st harmonic beats
                    fSynAdj - fOrb, ... % 1st harmonic beats
                    fSynAdj + fIoTA, ...
                    2*fSynAdj - fIoTA, ... % 2nd harmonic beats
                    2*fSynAdj + fIoTA, ...
                    3*fSynAdj - fIoTA, ... % 3rd harmonic beats
                    3*fSynAdj + fIoTA, ...
                    fOrbBeat, ... % ~ 2nd harmonic beat between orbital and TA periods
                    fSynAdj - fOrbMid2]; % Beat between strong(est) oscillations
                    

        case 'EUROPA'
            Pnut0 = Pnut0Jup_deg;
            Pnut1 = Pnut1Jup_degday;
            PnutM = [0, 0, 0, -0.980, -0.054, -0.014, 0.008, 0];
            fOrb = 101.3747235 / 360 / 86400;
            fTA = 1 / 84.62776 / 3600; %3.281745587e-6;
            fSyn = fJup - fOrb;
            fNearBeats = 1 / 3600 ./ [9.916007150192408304, 9.918866231039942249, ...
                9.924589341739000758, 12.954663079272373594];
            fNearOrbs = 1 / 3600 ./ [85.1513423, 84.575556990340160723, ...
                84.783999521416234302, 85.098596922214568394, 84.471719596866705615, ...
                84.888606553513909603, 84.993472034115839620, 85.415537694012940, ...
                85.521709897525156, 84.368136863961226, 85.628146376422208164];
            fOrbAdj = (101.3747235 + sum(PnutM .* Pnut1 .* cosd(Pnut0 + Pnut1*etMid_day))) / 360 / 86400;

            f = [fSyn, ... % synodic period
                 fTA, ... % true anomaly period,
                 fOrb]; % true orbital period wrt inertial reference frame

            % harmonics
            f = [f, 2*fSyn, ...
                    3*fSyn, ...
                    4*fSyn, ...
                    2*fTA, ...
                    2*fOrb, ...
                    fSyn - fTA, ... % 1st harmonic beats
                    fSyn + fTA, ...
                    fSyn - fOrb, ...
                    fSyn + fOrb, ...
                    2*fSyn - fTA, ... % 2nd harmonic beats
                    2*fSyn + fTA, ...
                    2*fTA - fOrb, ...
                    2*fOrb - fTA, ...
                    fSyn - 2*fOrb, ...
                    fSyn - 2*fTA, ...
                    fTA + fOrb, ...
                    fNearOrbs, ...
                    fNearBeats];


        case 'GANYMEDE'
            Pnut0 = Pnut0Jup_deg;
            Pnut1 = Pnut1Jup_degday;
            PnutM = [0, 0, 0, 0.033, -0.389, -0.082, 0, 0];
            fOrb = 50.3176081 / 360 / 86400; % orbital period
            fOrbAdj = (50.3176081 + sum(PnutM .* Pnut1 .* cosd(Pnut0 + Pnut1*etMid_day))) / 360 / 86400;
            fSyn = fJup - fOrb;
            fEuropaTA = 3.281745587e-6;
            fOrbEuropa = 101.3747235 / 360 / 86400;
            fMystery = 1 / 3600 / 34.724522173543868;
            fMystery2 = 1 / 3600 / 150.271894222476305458;
            f = [fSyn, fOrb];

            f = [f, 2*fSyn, ...
                    3*fSyn, ...
                    4*fSyn, ...
                    5*fSyn, ...
                    6*fSyn, ...
                    fSyn - fOrb, ... % 1st harmonic beats
                    fSyn + fOrb, ...
                    2*fSyn - fOrb, ... % 2nd harmonic beats
                    2*fSyn + fOrb, ...
                    3*fSyn - fOrb, ... % 3rd harmonic beats
                    3*fSyn + fOrb, ...
                    fMystery, ...
                    fSyn - fMystery, ...
                    fSyn + fMystery, ...
                    fEuropaTA/2, ... % Half Europa's true anomaly oscillation
                    fOrbAdj, ...
                    fJup - fOrbEuropa, ...
                    2*fSyn - 5*fOrb, ...
                    fMystery2, ...
                    fSyn - fMystery2, ...
                    1 / 3600 / 4.573237865242181, ...
                    1 / 3600 / 6.208585640926228, ...
                    1 / 3600 / 3.188814749737327, ...
                    1 / 3600 / 3.906250548182006, ...
                    2*fMystery];


        case 'CALLISTO'
            Pnut0 = Pnut0Jup_deg;
            Pnut1 = Pnut1Jup_degday;
            PnutM = [0, 0, 0, 0, 0.061, -0.533, 0, -0.009];
            fOrbAdj = (21.5710715 + sum(PnutM .* Pnut1 .* cosd(Pnut0 + Pnut1*etMid_day))) / 360 / 86400;
            fOrb = 21.5710715 / 360 / 86400;
            fSheet1 = 1 / 3600 / 15.120553206749488;
            fSheet2 = 1 / 3600 / 7.6696338880659605;
            fSheet3 = 1 / 3600 / 3.8072842207111743;
            fSheet4 = 1 / 3600 / 3.059002331397538;
            fSyn = fJup - fOrb;

            f = [fSyn, fOrbAdj];

            % harmonics
            f = [f, 2*fSyn, ...
                    3*fSyn, ...
                    4*fSyn, ...
                    5*fSyn, ...
                    6*fSyn, ...
                    2*fOrb, ...
                    fSyn - fOrb, ... % 1st harmonic beats
                    fSyn + fOrb, ...
                    2*fSyn - fOrb, ... % 2nd harmonic beats
                    2*fSyn + fOrb, ...
                    fSheet1, ...
                    fSheet2, ...
                    fSheet3, ...
                    fSheet4 ...
                    ];
            a=0;


        %% Saturn moons         
        case 'MIMAS'
            f = [1.23134297475e-5, ... % orbital period 1 
                 1.22487406465e-5, ... % orbital period 2
                 6.33328477817e-8]; % intermoon orbital resonance?

            % harmonics
            f = [f, 2*f(1), ...  
                    3*f(1), ...  
                    2*f(2), ...  
                    3*f(2)];  


        case 'ENCELADUS'
            fOrb = 262.7318996 / 360 / 86400; % orbital period
            fSyn = fSat - fOrb;
            fTA = 1 / 3600 / 32.927200612354675; % true anomaly (?) period
            f = [fOrb, fTA];

            f = [f, 2*fTA, ...
                    9.26559563e-6, ... % inter-moon positional resonances
                    4.62963117e-6, ... 
                    2.31798223e-6];

            f = [f, 1 / 3600 / 32.917867508555119116, ...
                    1 / 3600 / 32.941210199582918960, ...
                    1 / 3600 / 32.838749203627905615, ...
                    1 / 3600 / 32.936539015933050223, ...
                    1 / 3600 / 32.913202935703473884, ...
                    1 / 3600 / 32.834107031435593171, ...
                    1 / 3600 / 32.866629845374262686];


       case 'DIONE'
            f = [4.224300947e-6, ... % orbital period
                 4.2306342318e-6]; % orbital period


        %% Uranus moons   
        case 'MIRANDA'
            fSyn = 1/(35.057143696*3600);
            f = [fSyn, ... % synodic period
                 1/(33.9179560479*3600)]; % orbital period

            % harmonics
            f = [f, 2*f(1), ... % 1.5847e-05
                    3*f(1), ... % 2.3771e-05
                    4*f(1), ... % 3.1694e-05
                    2*f(2), ... % 1.6381e-5
                    f(2) - f(1), ... % 1st harmonic beats % 2.6613e-7  % this is large!
                    f(1) + f(2)]; % 1.6113e-5

            % 2nd harmonic beats
            f = [f, f(3) - f(2), ... % 7.6574e-06
                    f(3) + f(2), ... % 2.4305e-5
                    f(6) + f(1), ... % 2.3405e-5
                    f(4) - f(2), ... % 3rd harmonic beats % 1.5581e-05
                    f(4) + f(2), ... % 3.1960e-05
                    f(6) + f(3)]; % Double harmonic beats % 3.2228e-5        


        case 'ARIEL'
            fSyn = 1.15202411714e-5; % synodic period
            f(2) = 4.59162993363e-6; % orbital period

            % harmonics
            f = [f, 2*fSyn, ... % 6.028
                    3*fSyn, ... % 12.05
                    fSyn - f(2), ... % 1st harmonic beats % 40.09
                    fSyn + f(2)]; % 8.037


        case 'UMBRIEL'
            fSyn = 1.33189372223e-5;
            f = [fSyn, ... % synodic period
                 2.79273148465e-6]; % orbital period

            % harmonics
            f = [f, 2*fSyn, ...
                    3*fSyn, ...
                    4*fSyn, ...
                    fSyn - f(2), ... % 1st harmonic beats
                    fSyn + f(2)];


        case 'TITANIA'
            fSyn = 1.4782626327e-5; % synodic period
            f(2) = 1.32904237982e-6; % orbital period

            % harmonics
            f = [f, 2*fSyn, ...
                    3*fSyn, ...
                    4*fSyn, ...
                    fSyn - f(2), ... % 1st harmonic beats
                    fSyn + f(2)];


        case 'OBERON'
            fSyn = 1.52530978251e-5;
            f = [fSyn, ... % synodic period
                 8.60154960956e-7]; % orbital period

            % harmonics
            f = [f, 2*fSyn, ...
                    3*fSyn, ...
                    4*fSyn, ...
                    fSyn - f(2), ... % 1st harmonic beats
                    fSyn + f(2)];


        %% Neptune moons
        case 'TRITON'
            fSyn = 1.92125059641e-5;
            f = [fSyn, ... % synodic period
                 1.96941653582e-6]; % orbital period

            % harmonics
            f = [f, 2*fSyn, ...
                    3*fSyn, ...
                    2*f(2), ...
                    fSyn - f(2), ... % 1st harmonic beats
                    fSyn + f(2), ...
                    fSyn - 2*f(2), ...
                    fSyn + 2*f(2), ...
                    fSyn - 3*f(2), ...
                    fSyn + 3*f(2)];

            % 2nd harmonic beats
            f = [f, f(3) - f(2), ...
                    f(3) + f(2), ...
                    f(3) - 2*f(2), ...
                    f(3) + 2*f(2)];

    end

    f = sort(f, 'descend');
end

function Tout = T(f)
    Tout = 1 / 3600 ./ f;
end
