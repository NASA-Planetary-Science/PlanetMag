function [Brcs, Bpcs, Bzcs] = KS_TailMagNoTilt(Mode, rho, phi, z, localTime_rads)
% Calculate magnetotail field in current sheet coordinates.
%
% Parameters
% ----------
% Mode : int
%   Parameter selection for magnetotail model. Passing 7 or greater will include all coefficients.
% rho : double, 1xN
%   Distance from z axis in current sheet coordinates in planetary radii.
% phi : double, 1xN
%   Azimuthal angle in current sheet coordinates in radians.
% z : double, 1xN
%   Distance from current sheet normal plane in planetary radii.
% localTime_rads : double, 1xN
%   Local time in radians, i.e. east longitude in the JSO frame, which is the azimuthal angle
%   between the plane containing the Sun and the evaluation point.
%
% Returns
% -------
% Brcs, Bpcs, Bzcs : double, 1xN
%   Magnetic field contribution from magnetotail in current sheet cylindrical coordinates in nT.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Initialize outputs
    [Brcs, Bzcs] = deal(zeros(size(rho)));

    % Retrieve current sheet deformation model parameters
    [C, P, D0, D1, D, f, beta] = KS_coeffsCSdeform();
    
    % Get loop iteration for specified mode(s)
    if Mode < 7
        Modes = Mode:Mode;
    else
        Modes = 1:6;
    end
    
    drho = 0.05;
    dz = 0.05;
    
    for L=Modes
        zm = abs(z - dz);
        zp = abs(z + dz);
        % Adjust positions within current sheet half-thickness D
        zm(zm < D) = 0.5 * (D + zm(zm < D).^2/D);
        zp(zp < D) = 0.5 * (D + zp(zp < D).^2/D);
        
        % (Re-)initialize working variables
        [xlpp, xlpm] = deal(zeros(size(rho)));
    
        for i=1:6
            S1p = sqrt((beta(L,i) + zp).^2 + (rho + beta(L,i)).^2);
            S2p = sqrt((beta(L,i) + zp).^2 + (rho - beta(L,i)).^2);
            S1m = sqrt((beta(L,i) + zm).^2 + (rho + beta(L,i)).^2);
            S2m = sqrt((beta(L,i) + zm).^2 + (rho - beta(L,i)).^2);
            
            tp = 2*beta(L,i) ./ (S1p + S2p);
            tm = 2*beta(L,i) ./ (S1m + S2m);
            
            AAp = tp .* sqrt(1 - tp.^2) ./ (S1p .* S2p);
            AAm = tm .* sqrt(1 - tm.^2) ./ (S1m .* S2m);
            
            xlpp = xlpp + f(L,i) * AAp .* rho;
            xlpm = xlpm + f(L,i) * AAm .* rho;
        end
        
        Brcs = Brcs - (xlpp - xlpm) / (2*dz) * C(L);
    end
    
    for L=Modes
        rhom = rho - drho;
        rhop = rho + drho;
        % Simplifying from Khurana's code, zp and zm are used there with dz = 0 for both. We just
        % use zabs instead.
        zabs = abs(z);
        % Adjust positions within current sheet half-thickness D
        zabs(zabs < D) = 0.5 * (D + zabs(zabs < D).^2/D);
        
        % (Re-)initialize working variables
        [xlpp, xlpm] = deal(zeros(size(rho)));
        
        for i=1:6
            S1p = sqrt((beta(L,i) + zabs).^2 + (rhop + beta(L,i)).^2);
            S2p = sqrt((beta(L,i) + zabs).^2 + (rhop - beta(L,i)).^2);
            S1m = sqrt((beta(L,i) + zabs).^2 + (rhom + beta(L,i)).^2);
            S2m = sqrt((beta(L,i) + zabs).^2 + (rhom - beta(L,i)).^2);
            
            tp = 2*beta(L,i) ./ (S1p + S2p);
            tm = 2*beta(L,i) ./ (S1m + S2m);
            AAp = tp .* sqrt(1 - tp.^2) ./ (S1p .* S2p);
            AAm = tm .* sqrt(1 - tm.^2) ./ (S1m .* S2m);
            xlpp = xlpp + f(L,i) * AAp .* rhop;
            xlpm = xlpm + f(L,i) * AAm .* rhom;
        end
        
        Bzcs = Bzcs + (rhop.*xlpp - rhom.*xlpm) / (2*drho) ./ rho * C(L);
    end
    
    % Calculate Bphi from stretching parameter and radial field comp
    Bpcs = - P * rho .* Brcs;
end
