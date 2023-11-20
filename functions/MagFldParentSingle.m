function Bvec_nT = MagFldParentSingle(g, h, r_km, theta, phi, PlanetEqRadius, ...
    InternalFieldModel, ExternalFieldModel, magPhase_deg, Nmax, NmaxExt, SPHOUT)
% Evaluate the magnetic field of the desired planet at specified locations according to the
% specified internal and external field models.
%
% See MagFldPlanet for additional notes.
%
% Parameters
% ----------
% g : double, (Nmax)x(Nmax+1)
%   Internal field g coefficients in Schmidt normalization in G.
% h : double, (Nmax)x(Nmax+1)
%   Internal field h coefficients in Schmidt normalization in G.
% r_km : double
%   Radius of evaluation point in km from planet barycenter.
% theta : double
%   Planetocentric colatitude of evaluation points in radians.
% phi : double
%   Planetocentric east longitude of evaluation point in radians.
% PlanetEqRadius : double
%   Equatorial radius of planet in km associated with input ``g`` and ``h`` values.
% InternalFieldModel : char, 1xD
%   Code name for planetary intrinsic field model to evaluate. Refer to GetModelOpts for a listing
%   of valid options. Also accepts the special option ``'None'``, which skips evaluation of the
%   intrinsic field model.
% ExternalFieldModel : char, 1xE
%   Code name for current sheet field model to evaluate. Refer to GetModelOpts for a listing of
%   valid options. Also accepts the special option ``'None'``, which skips evaluation of the
%   current sheet field model.
% magPhase_deg : double, default=0
%   Arbitrary offset in degrees by which to rotate the magnetospheric field evaluation.
% Nmax : int, default=Inf
%   Maximum degree to which to limit intrinsic field models. Only has an effect if the value passed
%   is less than the lesser of the model maximum degree and the maximum implemented in spherical
%   harmonic calculations (currently 10).
% NmaxExt : int, default=Inf
%   Maximum degree to which to limit external field models. Only has an effect if the value passed
%   is less than the lesser of the model maximum degree and the maximum implemented in external
%   field calculations (currently 1).
% SPHOUT : bool, default=0
%   Whether to return vectors aligned to spherical coordinate axes (true) or cartesian (false).
%
% Returns
% -------
% Bvec_nT : double, 3x1
%   Magnetic field vector in standard coordinates (typically IAU or spherical coordinates if SPHOUT
%   is true) at evaluation points in nT. Output rows are x, y, z respectively for cartesian or r,
%   theta, phi for spherical.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~exist('SPHOUT', 'var'); SPHOUT = 0; end

    magPhase = deg2rad(magPhase_deg);
    npts = length(r_km);

    %% Adjust inputs and get dipole parameters
  
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

    r = r_km * 1e3;

    % Convert to cartesian coordinates, as needed for some calculations
    rxy = r .* sin(theta); % Projection onto xy plane (cylindrical coordinates)
    x = rxy .* cos(phi);
    y = rxy .* sin(phi);
    z = r .* cos(theta);

    % Inner Field
    if ~strcmp(InternalFieldModel,'None')

        dVr = 0; dVtheta = 0; dVphi = 0;
        dVrTemp = 0; dVthetaTemp = 0; dVphiTemp = 0;
        k = 0;
        for n = 1:Nmax
            k = k+1;
            for m = 0:k
                A = Rp_m*(Rp_m./r).^(n+1);
                dA = -(n+1)*(Rp_m./r).^(n+2);
                P = LegendreS(n,m,theta);
                dP = (1./r).*dLegendreS(n,m,theta);
                % m index of g and h are offset by 1 because Matlab cannot index 0
                Q = (g(n,m+1)*cos(m*phi)+h(n,m+1)*sin(m*phi));
                dQ = (1./(r.*sin(theta))) .* (-m*g(n,m+1)*sin(m*phi) + m*h(n,m+1)*cos(m*phi));

                ddVr = dA .* P .* Q;
                ddVtheta = A .* dP .* Q;
                ddVphi = A .* P .* dQ;

                save_order = '';
                % Save individual degree to show field line contribution
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

        % All but Khurana1997 and SphericalHarmonic current sheet models use the same basic design,
        % but with different params
        if ~(strcmp(ExternalFieldModel,'Khurana1997') || ...
            strcmp(ExternalFieldModel,'SphericalHarmonic'))

            if strcmp(ExternalFieldModel,'Connerney1981')
                opt = 2;
                if opt == 1
                    % Current sheet from 1981 publication, goes best with VIP4
                    Ri = 5; % Inner radius of current sheet in RJ
                    Ro = 50; % Outer radius of current sheet in RJ, initially 50
                    D = 2.5; % Half-thickness of current sheet in RJ
                    u0I0 = 0.0045; % Current constant in Gauss (Connerney 1982: u0I0/2 = 225 nT)
                    Theta0 = 9.6*pi/180; % Colatitude of sheet axis in rad
                    Phi0 = (360-202)*pi/180-magPhase;
                else
                    % Current sheet from JUNO workshop 2016
                    Ri = 5;
                    Ro = 56;
                    D = 3.1;
                    u0I0 = 0.0037;
                    Theta0 = 6.5*pi/180;
                    Phi0 = (360-206)*pi/180-magPhase;
                end

            elseif strcmp(ExternalFieldModel,'Connerney2020')
                % Current Sheet from Connerney 2020 JGR publication
                Ri = 7.8;
                Ro = 51.4;
                D = 3.6;
                Theta0 = 9.3*pi/180;
                Phi0 = (360-204.2)*pi/180-magPhase;
                %u0I0 = 0.003058; % Maximum (in paper, 152.9)
                % Average current constant in Gauss (Connerney 1982: u0I0/2 = 139.6 nT)
                %u0I0 = 0.002792;
                u0I0 = 0.002484; % Minimum (in paper, 124.2)

            elseif strcmp(ExternalFieldModel,'Cassini11')
                % Current Sheet from Dougherty et al. (2018) Table S2
                Ri = 6.5; % Inner radius of current sheet in RS
                Ro = 20; % Outer radius of current sheet in RS
                D = 2.5; % Half-thickness of current sheet in RS
                u0I0 = 0.0000400; % Current constant in Gauss for nominal model
                %u0I0 = 0.0000600; % Maximum from per-rev models
                %u0I0 = 0.0000299; % Minimum from per-rev models
                Theta0 = 0*pi/180; % Colatitude of sheet axis in rad
                Phi0 = 0*pi/180; % Longitude of sheet axis in rad
            end

            % Rotate to plasma sheet coordinates
            xm = x*cos(Phi0)*cos(Theta0) + y*sin(Phi0)*cos(Theta0) - z*sin(Theta0);
            ym = -x*sin(Phi0) + y*cos(Phi0);
            zm = x*cos(Phi0)*sin(Theta0) + y*sin(Phi0)*sin(Theta0) + z*cos(Theta0);

            % Convert to cylindrical coordinates divided by Rj
            rho = sqrt(xm.^2 + ym.^2) / Rp_m;
            psi = acos(xm/Rp_m./rho);
            psi(ym ~= 0) = psi.*sign(ym);
            zed = zm/Rp_m;

            % PlanetMag edit -- attenuate current sheet beyond certain distance.
            % No major effect to field lines within 50 RJ
            u0I0 = u0I0 * ones(1,npts);
            [eBrho, eBzed] = deal(zeros(1,npts));
            attentuation_distance = 50; % 50 RJ, way beyond Callisto orbit
            outer = find(rho > attentuation_distance);
            u0I0(outer) = u0I0(outer)./(rho(outer)/attentuation_distance).^2;   
            
            % Semi-infinite sheet with a = Ri
            a = Ri; % Inner radius
            
            % Region I: 0 < rho < 5RJ: within the start of the current sheet
            inner = find(rho < Ri);
            % a is held constant for p < inner radius, see last paragraph of Connerney 1981
            F1 = sqrt((zed(inner)-D).^2 + a^2);
            F2 = sqrt((zed(inner)+D).^2 + a^2);
            eBrho(inner) = (u0I0(inner)/2).*(rho(inner)/2).*(1./F1 - 1./F2);
            eBzed(inner) = (u0I0(inner)/2).*(2*D./sqrt((zed(inner).^2 + a^2)) ...
                - (rho(inner).^2/4).*((zed(inner)-D)./ F1.^3 - (zed(inner)+D)./F2.^3));
            % Below is as reported in the publication, but appears to result in a discontinuity in
            % field lines!
            % eBzed = (u0I0/2)*(2*D*(sqrt(zed^2 + a^2))^(-(1/a)/2) - (rho^2/4)*((zed-D)/ F1^3 ...
            %     - (zed+D)/F2^3));
            
            % Region II: rho > 5RJ AND Z > 2.5RJ (beyond start and above the current sheet)
            above = find(abs(zed) > D);
            F1 = sqrt((zed(above)-D).^2 + rho(above).^2);
            F2 = sqrt((zed(above)+D).^2 + rho(above).^2);
            eBrho(above) = (u0I0(above)/2).*((1./rho(above)).*(F1-F2+2*D*sign(zed(above))) ...
                - (a^2 .* rho(above)/4).*(1./F1.^3 - 1./F2.^3));
            eBzed(above) = (u0I0(above)/2).*(2*D./sqrt(zed(above).^2+rho(above).^2) ...
                - (a^2/4)*((zed(above)-D)./F1.^3 - (zed(above)+D)./F2.^3));

            % Region III: rho > 5RJ AND Z < 2.5RJ (within the current sheet)
            inside = find(~ ((rho < Ri) & (abs(zed) > D)) );
            F1 = sqrt((zed(inside)-D).^2 + rho(inside).^2);
            F2 = sqrt((zed(inside)+D).^2 + rho(inside).^2);
            eBrho(inside) = (u0I0(inside)/2).*((1./rho(inside)).*(F1-F2+2*zed(inside)) ...
                - (a^2*rho(inside)/4).*(1./F1.^3 -1./F2.^3));
            eBzed(inside) = (u0I0(inside)/2).*(2*D./sqrt(zed(inside).^2+rho(inside).^2) ...
                - (a^2/4)*((zed(inside)-D)./F1.^3 - (zed(inside)+D)./F2.^3)); % same as region II
            
            % Subtract from Semi-infinite sheet with a = Ro from Semi-infinite sheet with a = Ri
            % calculated above
            a = Ro; % Outer radius
            F1 = sqrt((zed-D).^2 + a^2);
            F2 = sqrt((zed+D).^2 + a^2);
            eBrho = eBrho - (u0I0/2).*(rho/2).*(1./F1 - 1./F2);
            eBzed = eBzed - (u0I0/2).*(2*D./sqrt((zed.^2 + a^2)) - (rho.^2/4).*((zed-D)./F1.^3 ...
                - (zed+D)./F2.^3));

            eBx =  eBrho.*cos(psi)*cos(Theta0)*cos(Phi0) + eBzed*sin(Theta0)*cos(Phi0) ...
                - eBrho.*sin(psi)*sin(Phi0);
            eBy =  eBrho.*cos(psi)*cos(Theta0)*sin(Phi0) + eBzed*sin(Theta0)*sin(Phi0) ...
                + eBrho.*sin(psi)*cos(Phi0);
            eBz = -eBrho.*cos(psi)*sin(Theta0) + eBzed*cos(Theta0);

        elseif strcmp(ExternalFieldModel,'Khurana1997') % Khurana plasma sheet model

            % Current sheet reference frame
            Theta0 = 9.6*pi/180; % Colatitude of sheet axis in rad
            % Longitude of sheet axis is 202 degrees SIII (1965) which is a left-handed coordinate
            % system. For IAU_JUPITER, 360-lambdaSIII for right-handed system
            Phi0 = (360-202)*pi/180-magPhase;

            thetacs = Theta0;

            % Constants: The best fit parameters obtained from Pioneer 10, Voyager 1, and Voyager 2
            % outbound data
            C1 = 80.3;
            C2 = 690.4;
            C3 = 101.3;
            C4 = -1.7;
            a1 = 2.49;
            a2 = 1.80;
            a3 = 2.64;
            r01 = 38.0;
            rho02 = 2.14;
            rho03 = 12.5;
            D1 = 2.01;
            D2 = 13.27;
            p = 6.26e-3;
            q = 0.35;
            x0 = -33.5;
            rho0 = 33.2;
            nu0 = 37.4;

            OmegaJ = 2*pi/9.92492; % Angular velocity of Jupiter in rad/hr

            % Convert to units of RJ
            xRj = x/Rp_m;
            yRj = y/Rp_m;
            zRj = z/Rp_m;
            rRj = r/Rp_m;

            % Convert to plasma sheet coordinates
            xm =  xRj*cos(Phi0)*cos(Theta0) + yRj*sin(Phi0)*cos(Theta0) - zRj*sin(Theta0);
            ym = -xRj*sin(Phi0) + yRj*cos(Phi0);
            zm =  xRj*cos(Phi0)*sin(Theta0) + yRj*sin(Phi0)*sin(Theta0) + zRj*cos(Theta0);
            rm =  rRj;

            % Convert to cylindrical coordinates
            rho = sqrt(xm.^2 + ym.^2);
            psi = acos(xm./rho);
            psi(ym ~= 0) = psi.*sign(ym);
            zed = zm;

            K0 = tan(thetacs);
            K1 = tanh(xm/x0);
            K2 = tanh(r01./rm);

            % Compute distance from plasma sheet
            delta = pi - (OmegaJ*rho0/nu0)*log(cosh(rho/rho0));
            % The distance between the current sheet and the jovigraphic equator at a cylindrical
            % radial distance of rho111 and the LH systemIII longitude of lambda
            Zcs = rho*K0.*((x0./xm).*K1.*cos(psi-delta)-cos(psi-pi));

            K3 = cosh((zed-Zcs)/D1);
            K4 = tanh((zed-Zcs)/D1);
            K5 = sech((zed-Zcs)/D2);
            K6 = tanh((zed-Zcs)/D2);

            % Compute partial derivatives
            dZcsdrho = K0*sech(xm/x0).^2 .* cos(psi-delta) ...
                - rho*K0.*(x0./xm).*K1.*sin(psi-delta).*(OmegaJ/nu0).*tanh(rho/rho0) ...
                - K0*cos(psi-pi);
            dZcsdpsi = -rho*K0.*((x0./xm).*K1.*sin(psi-delta)-sin(psi-pi));
            dfdrho = -C1 * K2.^a1 .* log(K3) ...
                + (C1*a1*r01.*rho.^2./rm.^3).*K2.^(a1-1).*sech(r01./rm).^2.*log(K3) ...
                + (C1*rho/D1).*K2.^a1 .* K4.*dZcsdrho + C2*rho.*tanh(rho02./rho).^a2 ...
                + C3*rho.*tanh(rho03./rho).^a3 + C4*rho;
            dfdpsi = (C1*rho/D1).*K2.^a1 .* K4.*dZcsdpsi;
            dfdzed = (C1*a1*r01.*rho.*zed./rm.^3) .* K2.^(a1-1) .* sech(r01./rm).^2 .* log(K3) ...
                - (C1*rho/D1).*K2.^a1 .* K4;
            dgdrho = p*(1+q*K6.^2) - (2*p*q*rho/D2).*K6.*K5.^2 .* dZcsdrho;
            dgdpsi = 1-(2*p*q*rho/D2).*K6.*K5.^2 .* dZcsdpsi;
            dgdzed = (2*p*q*rho/D2).*K6.*K5.^2;

            % Calculate field components from partial derivatives
            eBrho = (1./rho).*(dfdpsi.*dgdzed - dfdzed.*dgdpsi);
            eBpsi = dfdzed.*dgdrho - dfdrho.*dgdzed;
            eBzed = (1./rho).*(dfdrho.*dgdpsi - dfdpsi.*dgdrho);

            % Uncomment for zeroth order approximation
            % eBrho = - (1/rho)*dfdzed;
            % eBpsi = 0;
            % Bzed = (1/rho)*dfdrho;

            % Convert from cylindrical to cartesian coordinates
            eBxps = eBrho.*cos(psi) - eBpsi.*sin(psi);
            eByps = eBrho.*sin(psi) + eBpsi.*cos(psi);
            eBzps = eBzed;

            % Convert from plasma sheet coordinates to Jupiter coordinates
            eBx =  eBxps*cos(Theta0)*cos(Phi0) - eByps*sin(Phi0) + eBzps*sin(Theta0)*cos(Phi0);
            eBy =  eBxps*cos(Theta0)*sin(Phi0) + eByps*cos(Phi0) + eBzps*sin(Theta0)*sin(Phi0);
            eBz = -eBxps*sin(Theta0) + eBzps*cos(Theta0);

            % Convert from nT to Gauss for combining
            eBx = eBx * 1e-5;
            eBy = eBy * 1e-5;
            eBz = eBz * 1e-5;

        elseif strcmp(ExternalFieldModel,'SphericalHarmonic')
            
            dVr = 0; dVtheta = 0; dVphi = 0;
            k = 0;
            for n = 1:NmaxExt
                k = k+1;
                for m = 0:k

                    A = Rp_m*(r/Rp_m).^n;
                    dA = n*(r/Rp_m).^(n-1);

                    P = LegendreS(n,m,theta);
                    dP = (1./r).*dLegendreS(n,m,theta);
                    % m index of g and h are offset by 1 because Matlab cannot index 0
                    Q = (G(n,m+1)*cos(m*phi)+H(n,m+1)*sin(m*phi));
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
    
    % Output field in nT
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
