% Author: Corey J Cochrane and Marshall J Styczinski
% Inputs: r from body center (km), colatitude from angular momentum 
%   vector theta (rad), east longitude phi (rad) 
% Outputs: Bvec = [Bx (nT); By (nT); Bz (nT)], Magnetic Moment Mdip, and Dipole Offset Odip
% reference frame: IAU_PLANET (or US3 for Uranus), rotates with the planet

function [Bvec_nT, Mdip_nT, Odip_km] = MagFldParent(planet, r_km, theta, phi, InternalFieldModel, ...
                       ExternalFieldModel, magPhase_deg, SPHOUT, Nmaxin)
if ~exist('SPHOUT', 'var'); SPHOUT = 0; end
if ~exist('Nmaxin', 'var'); Nmaxin = Inf; end
magPhase = deg2rad(magPhase_deg); % J2000 phase offset for magnetic moment orientation (longitude)
npts = length(r_km);

[g, h, G, H, PlanetEqRadius, Nmax, NmaxExt] = GetGaussCoeffs(planet, InternalFieldModel);

    %% Adjust inputs and get dipole parameters
    
    % Limit Nmax if user desires
    Nmax = min(Nmax, Nmaxin);

    % Get planet radius in m
    Rp_m = PlanetEqRadius * 1000; % m
    
    % Adjust internal field coefficients if non-zero phase offset
    if magPhase ~= 0
        g_copy = g; h_copy = h;
        for n = 1:Nmax
            for m = 1:n % m = 0 element is axisymmetric and so does not change
                g(n,m+1) = g_copy(n,m+1) * cos(m*magPhase) - h_copy(n,m+1) * sin(m*magPhase);
                h(n,m+1) = h_copy(n,m+1) * cos(m*magPhase) + g_copy(n,m+1) * sin(m*magPhase);
            end
        end
    end

    % calculate dipole magnetic moment and dipole offset, ref: https://www.spenvis.oma.be/help/background/magfield/cd.html
    B0 = sqrt(g(1,1)^2+g(1,2)^2+h(1,2)^2); % reference field B0 = sqrt(g10^2+g11^2+h11^2), called reduced moment by Bartels
    M0 = 4*pi*B0*1e-15*Rp_m^3/(4e-7*pi); % dipole moment magnitude in SI units (Tm^3)
    thDip_rad = acos(g(1,1)/B0);  % dipole moment colatitude
    phiDip_rad = atan2(h(1,2), g(1,2)); % dipole moment east longitude
     % Dipole vector
    if SPHOUT
        Mdip_nT = [B0 * 1e5, thDip_rad, phiDip_rad];
    else
        Mdip_nT = [g(1,2), h(1,2), g(1,1)] * 1e5;
    end


    %% Offset dipole -- see Koochak and Fraser-Snith 2017
    % Note two typos in Koochak and Fraser-Snith (2017) Eq. 4: G11 and G20
    % should be g11 and g20, and the g20 inside square brackets should be
    % g21, see Fraser-Snith (1987).
    L0 = 2*g(1,1)*g(2,1) + sqrt(3)*(g(2,1)*g(2,2) + h(1,2)*h(2,2));                 % L0 = 2*g10*g20 + sqrt(3)*(g11*g21 + h11*h21)
    L1 =  -g(1,2)*g(2,1) + sqrt(3)*(g(1,1)*g(2,2) + g(1,2)*g(2,3) + h(1,2)*h(2,3)); % L1 =  -g11*g20 + sqrt(3)*(g10*g21 + g11*g22 + h11*h22)
    L2 =  -h(1,2)*g(2,1) + sqrt(3)*(g(1,1)*h(2,2) - h(1,2)*g(2,3) + g(1,2)*h(2,3)); % L2 =  -h11*g20 + sqrt(3)*(g10*h21 - h11*g22 + g11*h22)
    E = (L0*g(1,1) + L1*g(1,2) + L2*h(1,2)) / (4*B0^2);                             % E = (L0*g10 + L1*g11 + L2*h11) / (4*B0^2)

    % unitless offset parameters
    xi =   (L0 - g(1,1)*E) / (3*B0^2); % z-axis: xi =   (L0 - g10*E) / (3*B0^2)
    eta =  (L1 - g(1,2)*E) / (3*B0^2); % x-axis: eta =  (L1 - g11*E) / (3*B0^2)
    zeta = (L2 - h(1,2)*E) / (3*B0^2); % y-axis: zeta = (L2 - h11*E) / (3*B0^2)

    % Dipole offset in km
    if SPHOUT
        rO_km = sqrt(eta^2 + zeta^2 + xi^2);
        thO_rad = acos(xi / rO_km);
        phiO_rad = atan2(zeta, eta);
        Odip_km = PlanetEqRadius * [rO_km, thO_rad, phiO_rad];
    else
        Odip_km = PlanetEqRadius * [eta, zeta, xi];
    end
    
    r = r_km * 1e3;

    % Convert to cartesian coordinates, as needed for some calculations
    rxy = r .* sin(theta); % m, projection onto xy plane (cylindrical coordinates)
    x = rxy .* cos(phi); % m
    y = rxy .* sin(phi); % m
    z = r .* cos(theta); % m

    % Inner Field
    if ~strcmp(InternalFieldModel,'None')

        dVr = 0; dVtheta = 0; dVphi = 0; %radius, colatitude, longitude components
        dVrTemp = 0; dVthetaTemp = 0; dVphiTemp = 0;
        k = 0;
        for n = 1:Nmax  % degree, spherical harmonic index, n = 1 (dipole), n = 2 (quadrupole), n = 3 (octopole)
            k = k+1;
            for m = 0:k    % order
                A = Rp_m*(Rp_m./r).^(n+1);
                dA = -(n+1)*(Rp_m./r).^(n+2);
                P = LegendreS(n,m,theta);
                dP = (1./r).*dLegendreS(n,m,theta);
                Q = (g(n,m+1)*cos(m*phi)+h(n,m+1)*sin(m*phi));  % m index of g and h are offset by 1 because MATLAB cannot index 0
                dQ = (1./(r.*sin(theta))) .* (-m*g(n,m+1)*sin(m*phi) + m*h(n,m+1)*cos(m*phi));

                ddVr = dA .* P .* Q;
                ddVtheta = A .* dP .* Q;
                ddVphi = A .* P .* dQ;

                save_order = ''; % Undefined variable
                % save individual degree to show field line contribution
                if strcmp(save_order,'dipole')
                    if n == 1
                        dVrTemp = ddVr + dVrTemp;
                        dVthetaTemp = ddVtheta + dVthetaTemp;
                        dVphiTemp = ddVphi + dVphiTemp;
                    end
                elseif strcmp(save_order,'quadrupole')
                    if n == 2
                        dVrTemp =  ddVr + dVrTemp;
                        dVthetaTemp = ddVtheta + dVthetaTemp;
                        dVphiTemp = ddVphi + dVphiTemp;
                    end
                elseif strcmp(save_order,'octopole')
                    if n == 3
                        dVrTemp =  ddVr + dVrTemp;
                        dVthetaTemp = ddVtheta + dVthetaTemp;
                        dVphiTemp = ddVphi + dVphiTemp;
                    end
                end

                dVr = dVr + ddVr;
                dVtheta = dVtheta + ddVtheta;
                dVphi = dVphi + ddVphi;
            end

        end

        iBr = -dVr;
        iBth = -dVtheta;
        iBphi = -dVphi;

        [iBx, iBy, iBz] = Bsph2Bxyz(iBr, iBth, iBphi, theta, phi);

    else
        
        [iBx, iBy, iBz] = deal(zeros(size(r)));
    
    end

    % Outer Field
    if ~strcmp(ExternalFieldModel,'None')

        if ~(strcmp(ExternalFieldModel,'Khurana1997') || strcmp(ExternalFieldModel,'SphericalHarmonic'))  % All but Khurana1997 and SphericalHarmonic current sheet models use the same basic design, but with different params

            if strcmp(ExternalFieldModel,'Connerney1981')
                opt = 2;
                if opt == 1
                    % Current Sheet from 1981 publication, goes best with VIP4
                    Ri = 5;                       % Inner radius of current sheet in Rj
                    Ro = 50;                      % Outer radius of current sheet in Rj, initially 50
                    D = 2.5;                        % Half-thickness of current sheet in Rj
                    u0I0 = 0.0045;                % current constant in Gauss (Connerney 1982: u0I0/2 = 225nT)
                    Theta0 = 9.6*pi/180;          % Colatitude of sheet axis in rad
                    Phi0 = (360-202)*pi/180-magPhase;      % Longitude of sheet axis is 202 degrees SIII (1965) which is a left handed coordinate system ... for IAU_JUPITER, 360-lambdaSIII for right handed system
                else
                    % Current Sheet from JUNO workshop 2016
                    Ri = 5;                       % Inner radius of current sheet in Rj
                    Ro = 56;                      % Outer radius of current sheet in Rj, initially 50
                    D = 3.1;                        % Half-thickness of current sheet in Rj
                    u0I0 = 0.0037;                % current constant in Gauss (Connerney 1982: u0I0/2 = 225nT)
                    Theta0 = 6.5*pi/180;          % Colatitude of sheet axis in rad ... = (2/3) x 9.6 degrees
                    Phi0 = (360-206)*pi/180-magPhase;      % Longitude of sheet axis is 202 degrees SIII (1965) which is a left handed coordinate system ... for IAU_JUPITER, 360-lambdaSIII for right handed system
                end

            elseif strcmp(ExternalFieldModel,'Connerney2020')
                % Current Sheet from Connerney 2020 JGR publication
                Ri = 7.8;                       % Inner radius of current sheet in Rj
                Ro = 51.4;                      % Outer radius of current sheet in Rj, initially 50
                D = 3.6;                        % Half-thickness of current sheet in Rj
                Theta0 = 9.3*pi/180;          % Colatitude of sheet axis in rad ... = (2/3) x 9.6 degrees
                Phi0 = (360-204.2)*pi/180-magPhase;      % Longitude of sheet axis is 202 degrees SIII (1965) which is a left handed coordinate system ... for IAU_JUPITER, 360-lambdaSIII for right handed system
                %u0I0 = 0.003058;                % maximum (in paper, 152.9)
                %u0I0 = 0.002792;                % average current constant in Gauss (Connerney 1982: u0I0/2 = 139.6nT)
                u0I0 = 0.002484;                % minimum(in paper, 124.2)

            elseif strcmp(ExternalFieldModel,'Cassini11')
                % Current Sheet from Dougherty et al. (2018) Table S2
                Ri = 6.5;           % Inner radius of current sheet in Rs
                Ro = 20;            % Outer radius of current sheet in Rs, 
                D = 2.5;            % Half-thickness of current sheet in Rs
                u0I0 = 0.0000400;    % Current constant in Gauss for nominal model 
                %u0I0 = 0.0000299;        % Minimum from per-rev models
                %u0I0 = 0.0000600;        % Maximum from per-rev models
                Theta0 = 0*pi/180;  % Colatitude of sheet axis in rad
                Phi0 = 0*pi/180;    % Longitude of sheet axis in rad
            end

            % rotate to plasma sheet coordinates
            xm = x*cos(Phi0)*cos(Theta0) + y*sin(Phi0)*cos(Theta0) - z*sin(Theta0);
            ym = -x*sin(Phi0) + y*cos(Phi0);
            zm = x*cos(Phi0)*sin(Theta0) + y*sin(Phi0)*sin(Theta0) + z*cos(Theta0);

            % convert to cylindrical coordinates divided by Rj
            rho = sqrt(xm.^2 + ym.^2) / Rp_m;
            psi = acos(xm/Rp_m./rho);
            psi(ym ~= 0) = psi.*sign(ym);
            zed = zm/Rp_m;

            % PlanetMag edit -- attenuate current sheet beyond certain distance.
            % Visually no effect to field lines within 50 Rj
            u0I0 = u0I0 * ones(1,npts);
            [eBrho, eBzed] = deal(zeros(1,npts));
            attentuation_distance = 50; % 50 Rj way beyond Callisto orbit 
            outer = find(rho > attentuation_distance);
            u0I0(outer) = u0I0(outer)./(rho(outer)/attentuation_distance).^2;   
            
            % Semi-infinite sheet with a = Ri
            a = Ri; % inner radius
            
            % Region I: 0 < rho < 5Rj: within the start of the current sheet
            inner = find(rho < Ri);
            F1 = sqrt((zed(inner)-D).^2 + a^2); % a is held constant for p < inner radius, see last paragraph of Connerney 1981
            F2 = sqrt((zed(inner)+D).^2 + a^2);
            eBrho(inner) = (u0I0(inner)/2).*(rho(inner)/2).*(1./F1 - 1./F2);
            eBzed(inner) = (u0I0(inner)/2).*(2*D./sqrt((zed(inner).^2 + a^2)) - (rho(inner).^2/4).*((zed(inner)-D)./ F1.^3 - (zed(inner)+D)./F2.^3)); % - 2 * D / Ro);  % original code, appears to work properly
            %eBzed = (u0I0/2)*(2*D*(sqrt(zed^2 + a^2))^(-(1/a)/2) - (rho^2/4)*((zed-D)/ F1^3 - (zed+D)/F2^3));  % as reported in publication, appears to result with a discontinuity in field line!
            
            % Region II: rho > 5Rj AND Z > 2.5Rj (beyond of start and above the current sheet)
            above = find(abs(zed) > D);
            F1 = sqrt((zed(above)-D).^2 + rho(above).^2);
            F2 = sqrt((zed(above)+D).^2 + rho(above).^2);
            eBrho(above) = (u0I0(above)/2).*((1./rho(above)).*(F1-F2+2*D*sign(zed(above))) - (a^2 .* rho(above)/4).*(1./F1.^3 - 1./F2.^3));
            eBzed(above) = (u0I0(above)/2).*(2*D./sqrt(zed(above).^2+rho(above).^2) - (a^2/4)*((zed(above)-D)./F1.^3 - (zed(above)+D)./F2.^3));
            % Region III: rho > 5Rj AND Z < 2.5Rj (within the current sheet)
            
            inside = find(~ ((rho < Ri) & (abs(zed) > D)) );
            F1 = sqrt((zed(inside)-D).^2 + rho(inside).^2);
            F2 = sqrt((zed(inside)+D).^2 + rho(inside).^2);
            eBrho(inside) = (u0I0(inside)/2).*((1./rho(inside)).*(F1-F2+2*zed(inside)) - (a^2*rho(inside)/4).*(1./F1.^3 -1./F2.^3));
            eBzed(inside) = (u0I0(inside)/2).*(2*D./sqrt(zed(inside).^2+rho(inside).^2) - (a^2/4)*((zed(inside)-D)./F1.^3 - (zed(inside)+D)./F2.^3)); % same as region II
            
            % subtract from Semi-infinite sheet with a = Ro from Semi-infinite sheet with a = Ri calculated above
            a = Ro;  % outer radius
            F1 = sqrt((zed-D).^2 + a^2);
            F2 = sqrt((zed+D).^2 + a^2);
            eBrho = eBrho - (u0I0/2).*(rho/2).*(1./F1 - 1./F2);
            eBzed = eBzed - (u0I0/2).*(2*D./sqrt((zed.^2 + a^2)) - (rho.^2/4).*((zed-D)./F1.^3 - (zed+D)./F2.^3)); % - 2 * D / Ro);

            eBx =  eBrho.*cos(psi)*cos(Theta0)*cos(Phi0) + eBzed*sin(Theta0)*cos(Phi0) - eBrho.*sin(psi)*sin(Phi0);
            eBy =  eBrho.*cos(psi)*cos(Theta0)*sin(Phi0) + eBzed*sin(Theta0)*sin(Phi0) + eBrho.*sin(psi)*cos(Phi0);
            eBz = -eBrho.*cos(psi)*sin(Theta0) + eBzed*cos(Theta0);

        elseif strcmp(ExternalFieldModel,'Khurana1997')  % Khurana plasma sheet model

            % Current Sheet reference frame
            Theta0 = 9.6*pi/180;          % Colatitude of sheet axis in rad
            Phi0 = (360-202)*pi/180-magPhase;      % Longitude of sheet axis is 202 degrees SIII (1965) which is a left handed coordinate system ... for IAU_JUPITER, 360-lambdaSIII for right handed system

            thetacs = Theta0;

            % Constants: The Best Fit Parameters Obtained From Pioneer 10, Voyager I and Voyager 2 Outbound Data
            C1 = 80.3;
            C2 = 690.4;
            C3 = 101.3;
            C4 = -1.7;
            a1 = 2.49;
            a2 = 1.80;
            a3 = 2.64;
            r01 = 38.0;   % units Rj
            rho02 = 2.14; % units Rj
            rho03 = 12.5; % units Rj
            D1 = 2.01;    % units Rj
            D2 = 13.27;   % units Rj
            p = 6.26e-3;
            q = 0.35;
            x0 = -33.5;
            rho0 = 33.2;
            nu0 = 37.4;

            %OmegaJ = 1.76e-4*3600; % rad/hour, angular velocity of Jupiter, in rad/hr as nu0 has units of inverse hours
            OmegaJ = 2*pi/9.92492;  %  9.92492 hours

            % convert to units of Rj
            xRj = x/Rp_m;
            yRj = y/Rp_m;
            zRj = z/Rp_m;
            rRj = r/Rp_m;

            % convert to plasma sheet coordinates
            xm =  xRj*cos(Phi0)*cos(Theta0) + yRj*sin(Phi0)*cos(Theta0) - zRj*sin(Theta0);
            ym = -xRj*sin(Phi0) + yRj*cos(Phi0);
            zm =  xRj*cos(Phi0)*sin(Theta0) + yRj*sin(Phi0)*sin(Theta0) + zRj*cos(Theta0);
            rm =  rRj;

            %convert to cylindrical coordinates...
            rho = sqrt(xm.^2 + ym.^2);
            psi = acos(xm./rho);   % radians
            psi(ym ~= 0) = psi.*sign(ym);
            zed = zm;

            K0 = tan(thetacs);
            K1 = tanh(xm/x0);
            K2 = tanh(r01./rm);

            % compute distance from plasma sheet
            delta = pi - (OmegaJ*rho0/nu0)*log(cosh(rho/rho0));
            Zcs = rho*K0.*((x0./xm).*K1.*cos(psi-delta)-cos(psi-pi)); %The distance between the current sheet and the Jovigraphic equator at a cylindrical radial distance of rho111 and the LH systemIII longitude of lambda

            K3 = cosh((zed-Zcs)/D1);
            K4 = tanh((zed-Zcs)/D1);
            K5 = sech((zed-Zcs)/D2);
            K6 = tanh((zed-Zcs)/D2);

            % compute partial derivatives
            dZcsdrho = K0*sech(xm/x0).^2 .* cos(psi-delta) - rho*K0.*(x0./xm).*K1.*sin(psi-delta).*(OmegaJ/nu0).*tanh(rho/rho0) - K0*cos(psi-pi);
            dZcsdpsi = -rho*K0.*((x0./xm).*K1.*sin(psi-delta)-sin(psi-pi));
            dfdrho = -C1 * K2.^a1 .* log(K3) + (C1*a1*r01.*rho.^2./rm.^3).*K2.^(a1-1).*sech(r01./rm).^2.*log(K3) + (C1*rho/D1).*K2.^a1 .* K4.*dZcsdrho + C2*rho.*tanh(rho02./rho).^a2 + C3*rho.*tanh(rho03./rho).^a3 + C4*rho;
            dfdpsi = (C1*rho/D1).*K2.^a1 .* K4.*dZcsdpsi;
            dfdzed = (C1*a1*r01.*rho.*zed./rm.^3) .* K2.^(a1-1) .* sech(r01./rm).^2 .* log(K3) - (C1*rho/D1).*K2.^a1 .* K4;
            dgdrho = p*(1+q*K6.^2) - (2*p*q*rho/D2).*K6.*K5.^2 .* dZcsdrho;
            dgdpsi = 1-(2*p*q*rho/D2).*K6.*K5.^2 .* dZcsdpsi;
            dgdzed = (2*p*q*rho/D2).*K6.*K5.^2;

            % calculate field components from partial derivatives
            eBrho = (1./rho).*(dfdpsi.*dgdzed - dfdzed.*dgdpsi);
            eBpsi = dfdzed.*dgdrho - dfdrho.*dgdzed;
            eBzed = (1./rho).*(dfdrho.*dgdpsi - dfdpsi.*dgdrho);

            % uncomment for zeroth order approximation ...
            %eBrho = - (1/rho)*dfdzed;
            %eBpsi = 0;
            %Bzed = (1/rho)*dfdrho;

            % convert from cylindrical to cartesian coordinates
            eBxps = eBrho.*cos(psi) - eBpsi.*sin(psi);
            eByps = eBrho.*sin(psi) + eBpsi.*cos(psi);
            eBzps = eBzed;

            % convert from plasma sheet coordinates to Jupiter coordinates
            eBx =  eBxps*cos(Theta0)*cos(Phi0) - eByps*sin(Phi0) + eBzps*sin(Theta0)*cos(Phi0);
            eBy =  eBxps*cos(Theta0)*sin(Phi0) + eByps*cos(Phi0) + eBzps*sin(Theta0)*sin(Phi0);
            eBz = -eBxps*sin(Theta0) + eBzps*cos(Theta0);

            % convert from nT to Gauss for combining
            eBx = eBx * 1e-5;
            eBy = eBy * 1e-5;
            eBz = eBz * 1e-5;
            
            
        elseif strcmp(ExternalFieldModel,'SphericalHarmonic')
            
            dVr = 0; dVtheta = 0; dVphi = 0; % radius, colatitude, longitude components
            k = 0;
            for n = 1:NmaxExt  % degree, spherical harmonic index
                k = k+1;
                for m = 0:k    % order

                    A = Rp_m*(r/Rp_m).^n;
                    dA = n*(r/Rp_m).^(n-1);

                    P = LegendreS(n,m,theta);
                    dP = (1./r).*dLegendreS(n,m,theta);
                    Q = (G(n,m+1)*cos(m*phi)+H(n,m+1)*sin(m*phi));  % m index of g and h are offset by 1 because MATLAB cannot index 0
                    dQ = (1./(r.*sin(theta))) .* (-m*G(n,m+1)*sin(m*phi) + m*H(n,m+1)*cos(m*phi));

                    ddVr = dA .* P .* Q;
                    ddVtheta = A .* dP .* Q;
                    ddVphi = A .* P .* dQ;

                    dVr = dVr + ddVr;
                    dVtheta = dVtheta + ddVtheta;
                    dVphi = dVphi + ddVphi;
                end
            end

            eBr = -dVr;
            eBth = -dVtheta;
            eBphi = -dVphi;

            [eBx, eBy, eBz] = Bsph2Bxyz(eBr, eBth, eBphi, theta, phi);
            
        else
            
            error(['ExternalFieldModel ' ExternalFieldModel ' not recognized.']);
        
        end

    else
        
        [eBx, eBy, eBz] = deal(zeros(size(r)));
       
    end
    
    % output field in nT
    Bx = (iBx + eBx) * 1e5;
    By = (iBy + eBy) * 1e5;
    Bz = (iBz + eBz) * 1e5;
    
    % Optionally convert field vectors to spherical for output
    if SPHOUT
        [Br, Bth, Bphi] = Bxyz2Bsph(Bx, By, Bz, theta, phi);
        Bvec_nT = [Br; Bth; Bphi];
    else
        Bvec_nT = [Bx; By; Bz];
    end

end
