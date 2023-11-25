function [BrVIP4, BthVIP4, BphiVIP4] = KS_VIP4noDipole(r_Rp, theta, phi, gVIP4, hVIP4, AS_CODED)
% Calculate magnetic field vectors at evaluation points for the VIP4 model, sans the dipole moment.
%
% Get magnetic field vectors in System III spherical coordinates for the non-dipole multipole
% moments of the VIP4 model. We omit the dipole moments because we have already accounted for them
% in evaluating the shielded dipole magnetosphere model.
%
% Parameters
% ----------
% r_Rp : double, 1xN
%   Distance from Jupiter center of mass for evaluation points in terms of planet radius.
% theta : double, 1xN
%   Colatitude of evaluation points in System III coordinates in radians.
% phi : double, 1xN
%   East longitude of evaluation points in System III coordinates in radians.
% gVIP4 : double, 4x5
%   Schmidt semi-normalized g Gauss coefficient of VIP4 model in gauss.
% hVIP4 : double, 4x5
%   Schmidt semi-normalized h Gauss coefficient of VIP4 model in gauss.
% AS_CODED : bool, default=0
%   Whether to match the original Fortran code (true) or with increased precision and corrected
%   parameters (false). See MagFldJupiterKS2005 for more details.
%
% Returns
% -------
% BrVIP4, BthVIP4, BphiVIP4 : double, 1xN
%   Evalauted magnetic field vector components in System III spherical coordinates.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~AS_CODED
        % Convert distances to RJ_VIP4 = 71323 to fit with the model coefficients
        r_Rp = r_Rp * 71492 / 71323;

        % Copied from MagFldParent.m
        dVr = 0; dVtheta = 0; dVphi = 0;
        for n = 2:3
            for m = 0:n
                A = (1./r_Rp).^(n+1);
                dA = -(n+1)*(1./r_Rp).^(n+2);
                P = LegendreS(n,m,theta);
                dP = (1./r_Rp).*dLegendreS(n,m,theta);
                % m index of g and h are offset by 1 because Matlab cannot index 0
                Q = gVIP4(n,m+1)*cos(m*phi) + hVIP4(n,m+1)*sin(m*phi);
                dQ = (1./(r_Rp.*sin(theta))) .* (-m*gVIP4(n,m+1)*sin(m*phi) ...
                    + m*hVIP4(n,m+1)*cos(m*phi));

                dVr = dVr + dA .* P .* Q;
                dVtheta = dVtheta + A .* dP .* Q;
                dVphi = dVphi + A .* P .* dQ;
            end
        end

        BrVIP4 = -dVr;
        BthVIP4 = -dVtheta;
        BphiVIP4 = -dVphi;

    else
        % There are definitely simpler ways to perform these calculations, but in the interest of
        % trying to best match K. Khurana's code, we keep the same overall organization here.
        % MJS note: I tried and was unable to refactor this block of code, which is why the
        % alternate block above uses the spherical harmonic calculation from MagFldParent.

        rec = ones(1,91);
        for N=1:13
            N2 = (2*N - 1) * (2*N - 3);
            for M=1:N
                MN = N*(N-1)/2 + M;
                rec(MN) = (N-M) * (N+M-2) / N2;
            end
        end

        %        g0n,      g1n,      g2n,      g3n,     g4n
        G = [0,    0,        0, ...
            -0.05100, -0.61900,  0.49700, ...
            -0.01600, -0.52000,  0.24400, -0.17600, ...
            -0.16800,  0.22200, -0.06100, -0.20200, 0.06600, ...
             zeros(1,76)] * 1e5;
        %        h0n,      h1n,      h2n,      h3n,     h4n
        H = [0,    0,        0, ...
                   0, -0.36100,  0.05300, ...
                   0, -0.08800,  0.40800, -0.31600, ...
                   0,  0.07600,  0.40400, -0.16600, 0.03900, ...
             zeros(1,76)] * 1e5;
        S = 1;
        for N=2:13
            MN = N*(N-1)/2 + 1;
            S = S * (2*N-3) / (N-1);
            G(MN) = G(MN)*S;
            H(MN) = H(MN)*S;
            P = S;
            for M=2:N
                AA = 1;
                if M==2
                    AA = 2;
                end
                P = P * sqrt(AA * (N-M+1)/(N+M-2));
                MNN = MN + M - 1;
                G(MNN) = G(MNN) * P;
                H(MNN) = H(MNN) * P;
            end
        end

        PP = 1 ./ r_Rp;
        P = PP;
        [A, B] = deal(zeros(5, length(r_Rp)));
        for N=1:5
            P = P .* PP;
            A(N,:) = P;
            B(N,:) = P * N;
        end

        cos_phi = cos(phi);
        sin_phi = sin(phi);
        cos_th = cos(theta);
        sin_th = sin(theta);

        % Initialize working variables
        P = ones(size(r_Rp));
        [D, BBR, BBT, BBF] = deal(zeros(size(r_Rp)));

        for M=1:5
            if M == 1
                X = zeros(size(r_Rp));
                Y = ones(size(r_Rp));
            else
                MM = M - 1;
                W = X;
                X = W.*cos_phi + Y.*sin_phi;
                Y = Y.*cos_phi - W.*sin_phi;
            end

            Q = P;
            Z = D;
            [BI, P2, D2] = deal(zeros(size(r_Rp)));

            for N=M:5
                MN = N*(N-1)/2 + M;
                W = G(MN)*Y + H(MN)*X;

                P2(abs(P2) < 1e-38) = 0;
                Q(abs(Q) < 1e-38) = 0;

                BBR = BBR + B(N,:).*W.*Q;
                BBT = BBT - A(N,:).*W.*Z;

                if M ~= 1
                    QQ = Q;
                    QQ(sin_th < 1e-5) = Z(sin_th < 1e-5);
                    BI = BI + A(N,:).*QQ .* (G(MN)*X - H(MN)*Y);
                end
                DP = cos_th.*Z - sin_th.*Q - rec(MN)*D2;
                PM = cos_th.*Q - rec(MN).*P2;
                D2 = Z;
                P2 = Q;
                Z = DP;
                Q = PM;
            end

            D = sin_th.*D + cos_th.*P;
            P = sin_th .* P;

            if M ~= 1
                BI = BI * MM;
                BBF = BBF + BI;
            end
        end

        BrVIP4 = BBR;
        BthVIP4 = BBT;
        BphiVIP4 = BBF;
        BphiVIP4(sin_th > 1e-5) = BBF(sin_th > 1e-5) ./ sin_th(sin_th > 1e-5);
        BphiVIP4((sin_th < 1e-5) & (cos_th < 0)) = -BBF((sin_th < 1e-5) & (cos_th < 0));
    end
end
