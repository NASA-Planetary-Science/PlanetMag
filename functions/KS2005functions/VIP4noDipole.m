function [BrVIP4, BthVIP4, BphiVIP4] = VIP4noDipole(r, theta, phi, gVIP4, hVIP4)
% Get magnetic field vectors in System III spherical coordinates for the
% non-dipole multipole moments of the VIP4 model. We omit the dipole
% moments because we have already accounted for them in evaluating the
% shielded dipole magnetosphere model.

    global trueToKK;
    if ~trueToKK
        % Convert distances to Rj_VIP4 to fit with that model
        conv = 71492/71323; r = r * conv;
        Rp_m = 1.0;

        % Copied from MagFldJupiter.m
        dVr = 0; dVtheta = 0; dVphi = 0; %radius, colatitude, longitude components
        for n = 2:3  % degree, spherical harmonic index, n = 1 (dipole), n = 2 (quadrupole), n = 3 (octopole)
            for m = 0:n    % order
                A = Rp_m*(Rp_m./r).^(n+1);
                dA = -(n+1)*(Rp_m./r).^(n+2);
                P = LegendreS(n,m,theta);
                dP = (1./r).*dLegendreS(n,m,theta);
                Q = gVIP4(n,m+1)*cos(m*phi) + hVIP4(n,m+1)*sin(m*phi);  % m index of g and h are offset by 1 because MATLAB cannot index 0
                dQ = (1./(r.*sin(theta))) .* (-m*gVIP4(n,m+1)*sin(m*phi) + m*hVIP4(n,m+1)*cos(m*phi));

                dVr = dVr + dA .* P .* Q;
                dVtheta = dVtheta + A .* dP .* Q;
                dVphi = dVphi + A .* P .* dQ;
            end
        end

        BrVIP4 = -dVr;
        BthVIP4 = -dVtheta;
        BphiVIP4 = -dVphi;

    else
        % There are definitely simpler ways to perform these calculations,
        % but in the interest of trying to best match Khurana's code, we
        % keep the same overall organization here.

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

        PP = 1 ./ r;
        P = PP;
        [A, B] = deal(zeros(5, length(r)));
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
        P = ones(size(r));
        [D, BBR, BBT, BBF] = deal(zeros(size(r)));

        for M=1:5
            if M == 1
                X = zeros(size(r));
                Y = ones(size(r));
            else
                MM = M - 1;
                W = X;
                X = W.*cos_phi + Y.*sin_phi;
                Y = Y.*cos_phi - W.*sin_phi;
            end

            Q = P;
            Z = D;
            [BI, P2, D2] = deal(zeros(size(r)));

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
