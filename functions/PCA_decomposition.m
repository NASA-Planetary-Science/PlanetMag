function BD = PCA_decomposition(ets,body,Bx,By,Bz, descrip, PLOT_DIAGNOSTIC)
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
% input series (where Bi is one of Bx, By, or Bz). This calculation is made
% for each B component to get the excitation moments for each.

if ~exist('PLOT_DIAGNOSTIC', 'var')
    PLOT_DIAGNOSTIC = 1;
end
ets=ets(:); % Ensures ets is a row vector
npts = length(ets);
ets_ca = ets(round(length(ets)/2));
etMid_day = ets_ca / 86400;
%ets_ca = 1228738876.0; % for Triton flyby at AOL 50  degrees

fJup = 870.536 / 360 / 86400;
fSat = 810.7939024 / 360 / 86400;
fUra = abs(-501.1600928) / 360 / 86400;
fNep = 536.3128492 / 360 / 86400;
Pnut0Jup_deg = [73.32, 24.62, 283.90, 355.80, 119.90, 229.80, 352.25, 113.35];
Pnut1Jup_degday = [91472.9, 45137.2, 4850.7, 1191.3, 262.1, 64.3, 2382.6, 6070.0] ...
    / 36525;

switch body
    
    %% Jupiter moons       
    case 'EUROPA'
        parentName = 'Jupiter';
        % frequencies obtained from FFT : 5 year data computed from FFT of simulated data
        %fSyn = 2.47288637e-5; % synodic period
        %f(2) = 3.281745587e-6; % orbital period 1 .. stronger, higher freq
        %f(3) = 3.260134845e-6; % orbital period 2 
        
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
        %fOrbAdj = (101.3747235 + sum(PnutM .* Pnut1 .* cosd(Pnut0 + Pnut1*etMid_day))) / 360 / 86400;
        
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
        parentName = 'Jupiter';
        
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
        parentName = 'Jupiter';
        
        %f(1) = 2.729414e-5; % synodic period
        %f(2) = 6.939704e-7; % orbital period
        
        % good for callisto C09
        fSyn = 2.729408017371862e-05;
        f = [fSyn, ... % synodic period
             6.934070721821746e-07]; % orbital period
        
        %f(1) = 2.729713120183811e-05;
        %f(2) = 6.957100164384470e-07;
        
        % harmonics
        f = [f, 2*fSyn, ...
                3*fSyn, ...
                4*fSyn, ...
                5*fSyn, ...
                6*fSyn, ...
                7*fSyn, ...
                8*fSyn, ...
                9*fSyn, ...
                10*fSyn, ...
                2*f(2), ...
                fSyn - f(2), ... % 1st harmonic beats
                fSyn + f(2)];
        % 2nd harmonic beats
        f = [f, f(3) - f(2), ...
                f(3) + f(2), ...
                f(4) - f(2), ... % 3rd harmonic beats
                f(4) + f(2), ... % 4th harmonic beats
                f(5) - f(2), ...
                f(5) + f(2), ...
                fSyn - 2*f(2), ... % 1st harmonic double beats
                fSyn + 2*f(2)];
        
    
    %% Saturn moons         
    case 'MIMAS'
        parentName = 'Saturn';
        
        f = [1.23134297475e-5, ... % orbital period 1 
             1.22487406465e-5, ... % orbital period 2
             6.33328477817e-8]; % intermoon orbital resonance?
        
        % harmonics
        f = [f, 2*f(1), ...  
                3*f(1), ...  
                2*f(2), ...  
                3*f(2)];  
            
        
    case 'ENCELADUS'
        parentName = 'Saturn';
        
        fOrb = 262.7318996 / 360 / 86400; % orbital period
        fSyn = fSat - fOrb;
        fTA = 1 / 3600 / 32.927200612354675; % true anomaly (?) period
        f = [fOrb, fTA];
        
        f = [f, 2*fOrb, ...
                9.26559563e-6, ... % inter-moon positional resonances
                6.94761340e-6, ... 
                4.62963117e-6, ... 
                2.31798223e-6];
            
        
   case 'DIONE'
        parentName = 'Saturn';
        
        f = [4.224300947e-6, ... % orbital period
             4.2306342318e-6]; % orbital period
         
         
    %% Uranus moons   
    case 'MIRANDA'
        parentName = 'Saturn';
        
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
        parentName = 'Uranus';
        
        fSyn = 1.15202411714e-5; % synodic period
        f(2) = 4.59162993363e-6; % orbital period
        
        % harmonics
        f = [f, 2*fSyn, ... % 6.028
                3*fSyn, ... % 12.05
                fSyn - f(2), ... % 1st harmonic beats % 40.09
                fSyn + f(2)]; % 8.037
          
            
    case 'UMBRIEL'
        parentName = 'Uranus';
        
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
        parentName = 'Uranus';
        
        fSyn = 1.4782626327e-5; % synodic period
        f(2) = 1.32904237982e-6; % orbital period
        
        % harmonics
        f = [f, 2*fSyn, ...
                3*fSyn, ...
                4*fSyn, ...
                fSyn - f(2), ... % 1st harmonic beats
                fSyn + f(2)];
        
            
    case 'OBERON'
        parentName = 'Uranus';
        
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
        parentName = 'Neptune';
        
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
BdxEst = Bx*X*W;
BdyEst = By*X*W;
BdzEst = Bz*X*W;
Bxest = X*BdxEst';
Byest = X*BdyEst';
Bzest = X*BdzEst';

% Split up the excitation amplitudes into real (i for "in-phase") and
% imaginary (q for "quadrature") parts
Bdxi = BdxEst(1:2:end-2);
Bdxq = BdxEst(2:2:end-1);
BdxO = BdxEst(end);
Bdyi = BdyEst(1:2:end-2);
Bdyq = BdyEst(2:2:end-1);
BdyO = BdyEst(end);
Bdzi = BdzEst(1:2:end-2);
Bdzq = BdzEst(2:2:end-1);
BdzO = BdzEst(end);

% Fill output struct
BD.Bdxi = Bdxi;
BD.Bdxq = Bdxq;
BD.BdxO = BdxO;
BD.Bdyi = Bdyi;
BD.Bdyq = Bdyq;
BD.BdyO = BdyO;
BD.Bdzi = Bdzi;
BD.Bdzq = Bdzq;
BD.BdzO = BdzO;

BD.f = f;

disp(' ')
disp('Driving Field Waves...')
for i = 1:nFreqs
    disp([num2str(i), '. T = ', num2str(1/f(i)/3600,'%0.2f'), ' hr']) 
    Bdx = sqrt(Bdxi(i)^2+Bdxq(i)^2); Bdy = sqrt(Bdyi(i)^2+Bdyq(i)^2);Bdz = sqrt(Bdzi(i)^2+Bdzq(i)^2);
    %disp(['  2*pi*f*t = ', num2str(2*pi*f(i)*ets_ca*180/pi,'%0.2f'), ' deg'])
    disp(['  Bx = ', num2str(Bdx,'%0.2f'), ' nT, By = ', num2str(Bdy,'%0.2f'), ' nT, Bz = ', num2str(Bdz,'%0.2f'), ' nT'])
    disp(['  Px = ', num2str(atan2(Bdxq(i),Bdxi(i))*180/pi,'%0.2f'), ' deg, Py = ', num2str(atan2(Bdyq(i),Bdyi(i))*180/pi,'%0.2f'), ' deg, Pz = ', num2str(atan2(Bdzq(i),Bdzi(i))*180/pi,'%0.2f'), ' deg'])
end

disp(' ')

assignin('base','BdxEst',BdxEst)
assignin('base','BdyEst',BdyEst)
assignin('base','BdzEst',BdzEst)

%FieldCoeff.BdxEst = BdxEst;
%FieldCoeff.BdyEst = BdyEst;
%FieldCoeff.BdzEst = BdzEst;
%save 'FieldCoeff.mat'  FieldCoeff

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
        BxSyn = BdxO + Bdxi(f==fSyn)*cosSyn + Bdxq(f==fSyn)*sinSyn;
        BySyn = BdyO + Bdyi(f==fSyn)*cosSyn + Bdyq(f==fSyn)*sinSyn;
        BzSyn = BdzO + Bdzi(f==fSyn)*cosSyn + Bdzq(f==fSyn)*sinSyn;
    end

    figure; hold on; box on; grid on;
    set(gcf,'Name', [descrip ' model vs. reconstruction, first 200h']);
    title([bodyname ' ' descrip ' model vs. reconstruction, first 200h'])
    Bxin = plot(etsRel_h, Bx, 'b');
    Byin = plot(etsRel_h, By, 'r');
    Bzin = plot(etsRel_h, Bz, 'g');
    Bxout = plot(etsRel_h, Bxest, '--k');
    Byout = plot(etsRel_h, Byest, '--k');
    Bzout = plot(etsRel_h, Bzest, '--k');
    if ~strcmp(parentName,'Saturn')
        BxsynOnly = plot(etsRel_h, BxSyn, 'm');
        BysynOnly = plot(etsRel_h, BySyn, 'm');
        BzsynOnly = plot(etsRel_h, BzSyn, 'm');
        legend([Bxin,Byin,Bzin, Bxout, BxsynOnly], 'Bx', 'By', 'Bz', 'Reproduced', 'Synodic only')
    else
        legend([Bxin,Byin,Bzin, Bzout], 'Bx', 'By', 'Bz', 'Reproduced')
    end
    xlabel('Time (hr)')
    ylabel('Magnetic Field (nT)')
    set(gca,'fontsize',16)
    xlim([0 200])
    
    figure; hold on; box on; grid on;
    set(gcf,'Name', [descrip ' model vs. reconstruction diff, first 200h']);
    title([bodyname ' ' descrip ' model vs. reconstruction diff, first 200h'])
    Bxdiff = plot(etsRel_h, Bx-Bxest', 'b');
    Bydiff = plot(etsRel_h, By-Byest', 'r');
    Bzdiff = plot(etsRel_h, Bz-Bzest', 'g');
    if ~strcmp(parentName,'Saturn')
        Bxsyndiff = plot(etsRel_h, Bx-BxSyn', 'm');
        legend([Bxdiff,Bydiff,Bzdiff, Bxsyndiff], 'Bx','By','Bz', 'Bx synodic only')
    else
        legend([Bxdiff,Bydiff,Bzdiff], 'Bx','By','Bz')
    end
    xlabel('Time (hr)')
    ylabel('Magnetic Field Error (nT)')
    set(gca,'fontsize',16)
    xlim([0 200])
    
    figure; box on;
    set(gcf,'Name', [bodyname ' ' descrip ' hodogram']);
    if strcmp(parentName,'Saturn')
        plot(Bx, Bz, 'b')
        xlabel(['B_x IAU (\approx B_y ' bodyname(1) '\phi\Omega, nT)'])
        ylabel(['B_z IAU (\approx B_z ' bodyname(1) '\phi\Omega, nT)'])
        title([bodyname ' spin-parent plane hodogram, ' descrip])
    else
        plot(By, Bx, 'b')
        ylabel(['B_x IAU (\approx B_y ' bodyname(1) '\phi\Omega, nT)'])
        xlabel(['B_y IAU (\approx B_x ' bodyname(1) '\phi\Omega, nT)'])
        title([bodyname ' equatorial plane hodogram, ' descrip])
    end
    ylim([-max(abs(ylim())), max(abs(ylim()))])
    xlim(ylim())
    grid on;
    set(gca,'fontsize',16)
end

end

