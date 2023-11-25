function [Bperpr, Bperpf, Bperpx] = KS_BMPperp(rho, phi, x, Nmodes)
% Get magnetic field components perpendicular to the magnetopause boundary.
%
% The coordinates passed to this function are essentially cylindrical coordinates with the
% Jupiter--Sun--Orbital (JSO) x axis as the z axis of the coordinates and :math:`phi = 0` measured
% from the JSO y axis. See KS_S3CtoJSO for a definition of the JSO frame.
%
% Parameters
% ----------
% rho : double, 1xN
%   Distance from JSO x axis in planetary radii.
% phi : double, 1xN
%   Azimuthal angle between JSO y axis and evaluation point in radians.
% x : double, 1xN
%   Distance from JSO yz plane in planetary radii.
% Nmodes : int
%   "Number of dipole modes" is how the parameter was labeled in the original Fortran code. Its
%   usage in KS_BMPperp suggests it's some method of indexing Bessel functions and spherical
%   harmonics in that function.
%
% Returns
% -------
% Brds, Bpds, Bzds : double, 1xN
%   Magnetic field contribution from shielded dipole aligned to cylindrical System III axes in nT.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Get model parameters
    [A, B] = KS_coeffsMpause();
      
    Modes1 = Nmodes;
    Modes2 = 2 * Modes1;
    Modes3 = 3 * Modes1;
    Modes4 = 4 * Modes1;
    Modes5 = 5 * Modes1;
      
    % Initialize
    [Bperpr, Bperpf, Bperpx] = deal(zeros(size(rho)));
    
    rhom = rho - 0.1;
	rhop = rho + 0.1;
	phim = phi - 0.001;
	phip = phi + 0.001;
	phi2 = 2*phi;
	phi2m = 2*phi - 0.001;
	phi2p = 2*phi + 0.001;
	phi3 = 3*phi;
	phi3m = 3*phi - 0.001;
	phi3p = 3*phi + 0.001;
	xm = x - 0.1;
	xp = x + 0.1;
	sinf = sin(phi);
	cosf = cos(phi);
	sin2f = sin(phi2);
	cos2f = cos(phi2);
	sin3f = sin(phi3);
	cos3f = cos(phi3);
    sinp = sin(phip);
    cosp = cos(phip);
    sinm = sin(phim);
    cosm = cos(phim);
    sin2p = sin(phi2p);
    cos2p = cos(phi2p);
    sin2m = sin(phi2m);
    cos2m = cos(phi2m);
    sin3p = sin(phi3p);
    cos3p = cos(phi3p);
    sin3m = sin(phi3m);
    cos3m = cos(phi3m);
    
    % Calculate phi terms
    KA = 1;
    for KB=1:Modes1
        expx = exp(x / B(KB));
        expp = exp(xp / B(KB));
        expm = exp(xm / B(KB));
        Uperpp = besselj(1, rhop / B(KB)) .* expx .* sinf;
        Uperpm = besselj(1, rhom / B(KB)) .* expx .* sinf;
        Bperpr = Bperpr - A(KA)/0.2 * (Uperpp - Uperpm);
        
        Uperpp = besselj(1, rho / B(KB)) .* expx .* sinp;
        Uperpm = besselj(1, rho / B(KB)) .* expx .* sinm;
        Bperpf = Bperpf - A(KA)./(0.002*rho) .* (Uperpp - Uperpm);
        
        Uperpp = besselj(1, rho / B(KB)) .* expp .* sinf;
        Uperpm = besselj(1, rho / B(KB)) .* expm .* sinf;
        Bperpx = Bperpx - A(KA)/0.2 * (Uperpp - Uperpm);
        
        KA = KA + 1;
    end
    
    % Next set of A(KA) values for phi terms
    for KB=1:Modes1
        expx = exp(x / B(KB));
        expp = exp(xp / B(KB));
        expm = exp(xm / B(KB));
        Uperpp = besselj(1, rhop / B(KB)) .* expx .* cosf;
        Uperpm = besselj(1, rhom / B(KB)) .* expx .* cosf;
        Bperpr = Bperpr - A(KA)/0.2 * (Uperpp - Uperpm);
        
        Uperpp = besselj(1, rho / B(KB)) .* expx .* cosp;
        Uperpm = besselj(1, rho / B(KB)) .* expx .* cosm;
        Bperpf = Bperpf - A(KA)./(0.002*rho) .* (Uperpp - Uperpm);
        
        Uperpp = besselj(1, rho / B(KB)) .* expp .* cosf;
        Uperpm = besselj(1, rho / B(KB)) .* expm .* cosf;
        Bperpx = Bperpx - A(KA)/0.2 * (Uperpp - Uperpm);
        
        KA = KA + 1;
    end
    
    % Now calculate 2*phi terms
    for KB=(Modes1+1):Modes2
        expx = exp(x / B(KB));
        expp = exp(xp / B(KB));
        expm = exp(xm / B(KB));
        Uperpp = besselj(2, rhop / B(KB)) .* expx .* sin2f;
        Uperpm = besselj(2, rhom / B(KB)) .* expx .* sin2f;
        Bperpr = Bperpr - A(KA)/0.2 * (Uperpp - Uperpm);
        
        Uperpp = besselj(2, rho / B(KB)) .* expx .* sin2p;
        Uperpm = besselj(2, rho / B(KB)) .* expx .* sin2m;
        Bperpf = Bperpf - A(KA)./(0.002*rho) .* (Uperpp - Uperpm);
        
        Uperpp = besselj(2, rho / B(KB)) .* expp .* sin2f;
        Uperpm = besselj(2, rho / B(KB)) .* expm .* sin2f;
        Bperpx = Bperpx - A(KA)/0.2 * (Uperpp - Uperpm);
        
        KA = KA + 1;
    end
    
    % Next set of A(KA) values for 2*phi terms
    for KB=(Modes1+1):Modes2
        expx = exp(x / B(KB));
        expp = exp(xp / B(KB));
        expm = exp(xm / B(KB));
        Uperpp = besselj(2, rhop / B(KB)) .* expx .* cos2f;
        Uperpm = besselj(2, rhom / B(KB)) .* expx .* cos2f;
        Bperpr = Bperpr - A(KA)/0.2 * (Uperpp - Uperpm);
        
        Uperpp = besselj(2, rho / B(KB)) .* expx .* cos2p;
        Uperpm = besselj(2, rho / B(KB)) .* expx .* cos2m;
        Bperpf = Bperpf - A(KA)./(0.002*rho) .* (Uperpp - Uperpm);
        
        Uperpp = besselj(2, rho / B(KB)) .* expp .* cos2f;
        Uperpm = besselj(2, rho / B(KB)) .* expm .* cos2f;
        Bperpx = Bperpx - A(KA)/0.2 * (Uperpp - Uperpm);
        
        KA = KA + 1;
    end
    
    % Now calculate 3*phi terms
    for KB=(Modes2+1):Modes3
        expx = exp(x / B(KB));
        expp = exp(xp / B(KB));
        expm = exp(xm / B(KB));
        Uperpp = besselj(3, rhop / B(KB)) .* expx .* sin3f;
        Uperpm = besselj(3, rhom / B(KB)) .* expx .* sin3f;
        Bperpr = Bperpr - A(KA)/0.2 * (Uperpp - Uperpm);
        
        Uperpp = besselj(3, rho / B(KB)) .* expx .* sin3p;
        Uperpm = besselj(3, rho / B(KB)) .* expx .* sin3m;
        Bperpf = Bperpf - A(KA)./(0.002*rho) .* (Uperpp - Uperpm);
        
        Uperpp = besselj(3, rho / B(KB)) .* expp .* sin3f;
        Uperpm = besselj(3, rho / B(KB)) .* expm .* sin3f;
        Bperpx = Bperpx - A(KA)/0.2 * (Uperpp - Uperpm);
        
        KA = KA + 1;
    end
    
    % Next set of A(KA) values for 3*phi terms
    for KB=(Modes2+1):Modes3
        expx = exp(x / B(KB));
        expp = exp(xp / B(KB));
        expm = exp(xm / B(KB));
        Uperpp = besselj(3, rhop / B(KB)) .* expx .* cos3f;
        Uperpm = besselj(3, rhom / B(KB)) .* expx .* cos3f;
        Bperpr = Bperpr - A(KA)/0.2 * (Uperpp - Uperpm);
        
        Uperpp = besselj(3, rho / B(KB)) .* expx .* cos3p;
        Uperpm = besselj(3, rho / B(KB)) .* expx .* cos3m;
        Bperpf = Bperpf - A(KA)./(0.002*rho) .* (Uperpp - Uperpm);
        
        Uperpp = besselj(3, rho / B(KB)) .* expp .* cos3f;
        Uperpm = besselj(3, rho / B(KB)) .* expm .* cos3f;
        Bperpx = Bperpx - A(KA)/0.2 * (Uperpp - Uperpm);
        
        KA = KA + 1;
    end
    
    % Now calculate terms associated with derivative term and sin(phi)
    for KB=(Modes3+1):Modes4
        expx = exp(x / B(KB));
        expp = exp(xp / B(KB));
        expm = exp(xm / B(KB));
        rhob = rho / B(KB);
        besj0 = besselj(0, rhob);
        besj1 = besselj(1, rhob);
        Term1 = rhop.*besselj(0, rhop/B(KB)) + (x - B(KB)).*besselj(1, rhop/B(KB));
        Term2 = rhom.*besselj(0, rhom/B(KB)) + (x - B(KB)).*besselj(1, rhom/B(KB));
        Term3 = expx .* (rho.*besj0 + (x - B(KB)).*besj1);
        Term4 = expp .* (rho.*besj0 + (xp - B(KB)).*besj1);
        Term5 = expm .* (rho.*besj0 + (xm - B(KB)).*besj1);
        
        Uperpp = Term1 .* expx .* sinf;
        Uperpm = Term2 .* expx .* sinf;
        Bperpr = Bperpr - A(KA)/0.2 * (Uperpp - Uperpm);
        Uperpp = Term3 .* sinp;
        Uperpm = Term3 .* sinm;
        Bperpf = Bperpf - A(KA)./(0.002*rho) .* (Uperpp - Uperpm);
        Uperpp = Term4 .* sinf;
        Uperpm = Term5 .* sinf;
        Bperpx = Bperpx - A(KA)/0.2 * (Uperpp - Uperpm);
        
        KA = KA + 1;
    end
    
    % Now the same for derivative term and cos(phi)
    for KB=(Modes3+1):Modes4
        expx = exp(x / B(KB));
        expp = exp(xp / B(KB));
        expm = exp(xm / B(KB));
        rhob = rho / B(KB);
        besj0 = besselj(0, rhob);
        besj1 = besselj(1, rhob);
        Term1 = rhop.*besselj(0, rhop/B(KB)) + (x - B(KB)).*besselj(1, rhop/B(KB));
        Term2 = rhom.*besselj(0, rhom/B(KB)) + (x - B(KB)).*besselj(1, rhom/B(KB));
        Term3 = expx .* (rho.*besj0 + (x - B(KB)).*besj1);
        Term4 = expp .* (rho.*besj0 + (xp - B(KB)).*besj1);
        Term5 = expm .* (rho.*besj0 + (xm - B(KB)).*besj1);
        
        Uperpp = Term1 .* expx .* cosf;
        Uperpm = Term2 .* expx .* cosf;
        Bperpr = Bperpr - A(KA)/0.2 * (Uperpp - Uperpm);
        Uperpp = Term3 .* cosp;
        Uperpm = Term3 .* cosm;
        Bperpf = Bperpf - A(KA)./(0.002*rho) .* (Uperpp - Uperpm);
        Uperpp = Term4 .* cosf;
        Uperpm = Term5 .* cosf;
        Bperpx = Bperpx - A(KA)/0.2 * (Uperpp - Uperpm);
        
        KA = KA + 1;
    end
    
    % Now calculate terms associated with 2nd derivative term and sin(phi), or something like that
    for KB=(Modes4+1):Modes5
        expx = exp(x / B(KB));
        expp = exp(xp / B(KB));
        expm = exp(xm / B(KB));
        rhob = rho / B(KB);
        besj2 = besselj(2, rhob);
        besj3 = besselj(3, rhob);
        Term1 = rhop.*besselj(2, rhop/B(KB)) + (x - 3*B(KB)).*besselj(3, rhop/B(KB));
        Term2 = rhom.*besselj(2, rhom/B(KB)) + (x - 3*B(KB)).*besselj(3, rhom/B(KB));
        Term3 = expx .* (rho.*besj2 + (x - 3*B(KB)).*besj3);
        Term4 = expp .* (rho.*besj2 + (xp - 3*B(KB)).*besj3);
        Term5 = expm .* (rho.*besj2 + (xm - 3*B(KB)).*besj3);
        
        Uperpp = Term1 .* expx .* sin3f;
        Uperpm = Term2 .* expx .* sin3f;
        Bperpr = Bperpr - A(KA)/0.2 * (Uperpp - Uperpm);
        Uperpp = Term3 .* sin3p;
        Uperpm = Term3 .* sin3m;
        Bperpf = Bperpf - A(KA)./(0.002*rho) .* (Uperpp - Uperpm);
        Uperpp = Term4 .* sin3f;
        Uperpm = Term5 .* sin3f;
        Bperpx = Bperpx - A(KA)/0.2 * (Uperpp - Uperpm);
        
        KA = KA + 1;
    end
    
    % Same for cos(phi)
    for KB=(Modes4+1):Modes5
        expx = exp(x / B(KB));
        expp = exp(xp / B(KB));
        expm = exp(xm / B(KB));
        rhob = rho / B(KB);
        besj2 = besselj(2, rhob);
        besj3 = besselj(3, rhob);
        Term1 = rhop.*besselj(2, rhop/B(KB)) + (x - 3*B(KB)).*besselj(3, rhop/B(KB));
        Term2 = rhom.*besselj(2, rhom/B(KB)) + (x - 3*B(KB)).*besselj(3, rhom/B(KB));
        Term3 = expx .* (rho.*besj2 + (x - 3*B(KB)).*besj3);
        Term4 = expp .* (rho.*besj2 + (xp - 3*B(KB)).*besj3);
        Term5 = expm .* (rho.*besj2 + (xm - 3*B(KB)).*besj3);
        
        Uperpp = Term1 .* expx .* cos3f;
        Uperpm = Term2 .* expx .* cos3f;
        Bperpr = Bperpr - A(KA)/0.2 * (Uperpp - Uperpm);
        Uperpp = Term3 .* cos3p;
        Uperpm = Term3 .* cos3m;
        Bperpf = Bperpf - A(KA)./(0.002*rho) .* (Uperpp - Uperpm);
        Uperpp = Term4 .* cos3f;
        Uperpm = Term5 .* cos3f;
        Bperpx = Bperpx - A(KA)/0.2 * (Uperpp - Uperpm);
        
        KA = KA + 1;
    end
    
    % Negate for output
    Bperpr = -Bperpr;
    Bperpf = -Bperpf;
    Bperpx = -Bperpx;
end
