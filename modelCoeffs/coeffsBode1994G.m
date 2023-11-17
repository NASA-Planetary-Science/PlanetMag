function Gnm = coeffsBode1994G(alpha)
% **Gnm**
%
% Calculate Jupiter magnetopause field coefficients as fit by Bode (1994).
%
% For use with Engle (1992) method as a function of dipole "precession angle" alpha for Jupiter.
% See https://apps.dtic.mil/sti/pdfs/ADA284857.pdf.
%
% Parameters
% ----------
% alpha : double, 1xN
%   Nod longitude for planetary dipole moment in degrees. This is the angle between the equatorial
%   projection of the magnetic dipole moment vector and the direction of the Sun.
%   alpha = 0 is when the the dipole nod longitude is at local noon.
%
% Returns
% -------
% Gnm : double, 10x11xN
%   Numerical coefficients that describe the fit by Bode (1994) for each dipole angle alpha in the
%   time series to be modeled.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    npts = length(alpha);
    Gnm = zeros(10, 11, npts);
    fac1 = 0.38/pi^2;
    fac3 = 0.32/9/pi^2;
    fac5 = 0.32/25/pi^2;
    fac7 = 0.4/49/pi^2;
    cos1 = cos(alpha);
    cos2 = cos(2*alpha);
    cos3 = cos(3*alpha);
    cos4 = cos(4*alpha);
    cos5 = cos(5*alpha);
    cos7 = cos(7*alpha);
    sin05 = sin(0.5*alpha);
    
    Gnm(1,1,:) = -(0.0403*cos1 - 0.98332*fac1*cos1 - 3.10310*fac3*cos3 - ...
        0.806*fac5*cos5 - 0.806*fac7*cos7 + 0.00702*cos2 + 0.00494*cos4 - ...
        0.0172*sin05 - 0.5823);
    Gnm(1,2,:) = -(6*(fac1*cos1 + fac3*cos3 + fac5*cos5 + fac7*cos7) - ...
        0.35*cos1 - 0.008*cos2 + 0.007*cos4 - 0.0005);
    Gnm(2,1,:) = -(2.475*(fac1*cos1 + fac3*cos3 + fac5*cos5 + ...
        fac7*cos7) - 0.97*cos1 - 0.0032*cos4 + 0.003*cos2 - fac5/3.2*cos5 + 0.001);
    Gnm(2,2,:) = -(0.0969*cos1 + 0.13357*sin05 -0.855*(fac1*cos1 + fac3*cos3 + ...
        fac5*cos5 + fac7*cos7) + 0.00551*cos1 - 0.00285*cos4 - 0.41);
    Gnm(2,3,:) = -(0.3*(fac1*cos1 + fac3*cos3 + fac5*cos5 + ...
        fac7*cos7) - 0.0593*cos1 - 0.0029*cos2 - 0.00158*cos4 + 0.0047);
    Gnm(3,1,:) = -(-0.03773*cos1 - 0.04921*sin05 + 0.287*(fac1*cos1 + fac3*cos3 + ...
        fac5*cos5 + fac7*cos7) + 0.0063*cos4 - 0.0077*cos2 - 0.051);
    Gnm(3,2,:) = -(3.1*(fac1*cos1 + fac3*cos3 + fac5*cos5 + ...
        fac7*cos7) - 0.1835*cos1 - 0.0001);
    Gnm(3,3,:) = -(0.05165*cos1 + 0.1*sin05 - 0.1325*(fac1*cos1 + fac3*cos3 + ...
        fac5*cos5 + fac7*cos7) + 0.0029*cos1 + 0.0013*cos4 - 0.1268);
    Gnm(3,4,:) = -(-0.3096*(fac1*cos1 + fac3*cos3 + fac5*cos5 + fac7*cos7) + 0.001);
    Gnm(4,1,:) = -(0.9625*(fac1*cos1 + fac3*cos3 + fac5*cos5 - ...
        fac7/2*cos7 + 0.6/16/pi^2*cos2 - 0.24/16/pi^2*cos4) - 0.00225);
    Gnm(4,2,:) = -(0.024012*cos1 + 0.046*sin05 - 0.00276*cos2 - 0.0092*cos4 - 0.16074);
    Gnm(4,3,:) = -(3.415*(fac1*cos1 + fac3*cos3 + fac5*cos5 + ...
        fac7*cos7) - 0.181*cos1 - 0.002*cos2 + 0.0025);
    Gnm(4,4,:) = -(0.8*(-0.005*cos2 + 0.00122*cos1 + 0.001*cos4) + 0.0003);
    Gnm(4,5,:) = -(0.65*(0.56*(fac1*cos1 + fac3*7/3.2*cos3 + fac5*cos5 + ...
        fac7*cos7) - 0.0261*cos1 - 0.003*cos2) - 0.001);
    Gnm(5,1,:) = -(0.0025*cos3 + 0.00084*cos1 - 0.00072*cos2 + 0.0001*cos4 - 0.00715);
    Gnm(5,2,:) = -(0.5*(1.6*(fac1*cos1 + fac3*95/32*cos3 + fac5*cos5 + ...
        fac7*cos7) - 0.0985*cos1 + 0.004*cos4 + 0.002*cos2) - 0.005);
    Gnm(5,3,:) = -(0.006*cos3 + 0.008*cos2 - 0.0034*cos1 - 0.002*cos4 - 0.001*cos5 - 0.0695);
    Gnm(5,4,:) = -(0.6*(-0.0303*cos1 + 0.03*cos3 - 0.006*cos4) + 0.0052);
    Gnm(5,5,:) = -(-0.00296*cos1 + 0.0033*cos2 + 0.001*cos4 - 0.0019);
    Gnm(5,6,:) = -(0.77*(-0.0022*cos1 + 0.0042*cos3 - 0.003*cos2 + 0.0006*cos4 + 0.0013));
    Gnm(6,1,:) = -(0.279248671/pi^2*cos1 - 0.0002868*cos2 + 0.02936985/9/pi^2*cos3 - ...
        0.00346*cos4 + 0.193610650/25/pi^2*cos5 + 0.193610650/49/pi^2*cos7 + 0.003016240);
    Gnm(6,2,:) = -(0.0082*cos1 - 0.0103*cos2 - 0.0092*cos4 - 0.0517);
    Gnm(6,3,:) = -(-0.028*cos1 - 0.0035*cos2 + 0.0072*cos3 + 0.0015*cos4 - 0.0018);
    Gnm(6,4,:) = -(0.003*cos1 + 0.017*cos2 + 0.0017*cos4 - 0.0225);
    Gnm(6,5,:) = -(0.72*(-0.0058*cos1 - 0.0051*cos2 + 0.015*cos3 - 0.007*cos4) + 0.001);
    Gnm(6,6,:) = -(0.485*(-0.0037*cos1 + 0.013*cos2 - 0.0075*cos3) - 0.0065);
    Gnm(6,7,:) = -(0.81*(-0.00235*cos1 - 0.0022*cos2 + 0.004*cos3 - 0.00032*cos4) + 0.00181);
    Gnm(7,1,:) = -(0.72*(-0.007*cos1 + 0.005*cos2 + 0.009*cos3 + 0.015*cos4) + 0.007);
    Gnm(7,2,:) = -(0.92*(0.0041*cos1 + 0.01*cos2 + 0.0036*cos4) + 0.0007);
    Gnm(7,3,:) = -(1.03*(0.0022*cos1 - 0.0075*cos2 + 0.0034*cos3 - 0.01*cos4) - 0.036);
    Gnm(7,4,:) = -(0.83*(-0.0204*cos1 - 0.001*cos2 + 0.017*cos3 - 0.0027*cos4) + 0.0032);
    Gnm(7,5,:) = -(0.955*(-0.0024*cos1 + 0.012*cos2 + 0.0012*cos3 + 0.0007*cos4) - 0.00135);
    Gnm(7,6,:) = -(0.815*(-0.0013*cos1 - 0.0041*cos2 + 0.006*cos3 - 0.0078*cos4) + 0.0013);
    Gnm(7,7,:) = -(0.94*(0.0042*cos1 + 0.0044*cos2 - 0.006*cos3 + 0.0025*cos4) - 0.0063);
    Gnm(7,8,:) = -(0.6*(-0.00335*cos1 - 0.001*cos2 + 0.005*cos3 - 0.0015*cos4) + 0.00035);
    Gnm(8,1,:) = -(0.92*(0.0065*cos1 -0.008*cos2 + 0.0032*cos3 - 0.0032*cos4) - 0.0008);
    Gnm(8,2,:) = -(0.9*(0.000013*cos1 + 0.00005*cos3 + 0.00008*cos4) - 0.0000237);
    Gnm(8,3,:) = -(0.6*(-0.00005*cos1 + 0.0001*cos2 - 0.00012*cos3 + 0.000075*cos4) - 0.000015);
    Gnm(8,4,:) = -(0.665*(0.005*cos1 - 0.0005*cos2 + 0.003*cos3 - 0.0035*cos4) - 0.019);
    Gnm(8,5,:) = -(-0.00207*cos1 - 0.00105*cos2 + 0.005*cos3 - 0.0006*cos4 + 0.00012);
    Gnm(8,6,:) = -(0.88*(-0.0019*cos1 + -0.005*cos2 - 0.0019*cos3) - 0.00085);
    Gnm(8,7,:) = -(0.9*(0.0022*cos1 - 0.004*cos2 + 0.0022*cos3 - 0.003*cos4) + 0.0005);
    Gnm(8,8,:) = -(0.8*(0.0036*cos1 + 0.001*cos2 - 0.005*cos3 + 0.004*cos4) - 0.00288);
    Gnm(8,9,:) = -(-0.0006*cos1 - 0.0009*cos2 + 0.0014*cos3 + 0.0005*cos4 - 0.00085);
    Gnm(9,1,:) = -(1.1*(-0.0022*cos1 - 0.002*cos2 + 0.0137*cos3 + 0.015*cos4) - 0.012);
    Gnm(9,2,:) = -(1.33*(-0.0001*cos1 + 0.00033*cos2 - 0.0005*cos3 + 0.0002*cos4) + 0.00003);
    Gnm(9,3,:) = -(1.05*(-0.00016*cos1 + 0.00024*cos2 + 0.00006*cos3 - 0.00035*cos4) - 0.00002);
    Gnm(9,4,:) = -(0.95*(-0.00285*cos1 - 0.0037*cos2 + 0.005*cos3 - 0.0003*cos4 + 0.0010575));
    Gnm(9,5,:) = -(0.75*(-0.0026*cos1 + 0.0058*cos2 + 0.005*cos3 - 0.0032*cos4) - 0.0048);
    Gnm(9,6,:) = -(0.91*(-0.0024*cos1 + 0.0009*cos2 + 0.004*cos3 - 0.0022*cos4) - 0.00008);
    Gnm(9,7,:) = -(0.7*(0.000515*cos1 + 0.00071*cos2 + 0.0006*cos3 - 0.00032*cos4) - 0.00158);
    Gnm(9,8,:) = -(1.05*(0.00154*cos1 - 0.0018*cos2 + 0.0003*cos3 - 0.0011*cos4) + 0.00085);
    Gnm(9,9,:) = -(1.64*(0.00061*cos1 + 0.00001*cos2 - 0.001*cos3 + 0.0008*cos4) - 0.00069);
    Gnm(9,10,:) = -(1.13*(0.001*cos1 - 0.00145*cos2 + 0.0003*cos3 + 0.00055*cos4) - 0.00077);
    Gnm(10,1,:) = -(0.7*(0.00019*cos1 - 0.0004*cos2 + 0.0005*cos3 - 0.00026*cos4) + 0.000002);
    Gnm(10,2,:) = -(0.55*(-0.00004*cos1 - 0.00024*cos2 + 0.0004*cos3 + 0.0005*cos4) - 0.00021);
    Gnm(10,3,:) = -(0.88*(-0.00004*cos1 - 0.00024*cos2 + 0.0004*cos3 + 0.0005*cos4) - 0.00021);
    Gnm(10,4,:) = -(0.00018*cos2 + 0.00083*cos3 - 0.0001*cos4 - 0.00157);
    Gnm(10,5,:) = -(1.16*(0.0008*cos1 - 0.00125*cos2 + 0.0006*cos3 - 0.00002*cos4) - 0.00007);
    Gnm(10,6,:) = -(-0.00185*cos1 + 0.00185*cos2 + 0.0001*cos3 - 0.0002*cos4 + 0.00034);
    Gnm(10,7,:) = -(0.00045*cos1 + 0.00013*cos2 + 0.0007*cos3 - 0.00063*cos4 - 0.00071);
    Gnm(10,8,:) = -(0.00022*cos1 + 0.00026*cos2 + 0.00058*cos3 - 0.00067*cos4 - 0.00039);
    Gnm(10,9,:) = -(0.0018*cos1 - 0.0014*cos2 - 0.0005*cos3 - 0.0003*cos4 + 0.0004);
    Gnm(10,10,:) = -(1.04*(-0.0000026*cos1 - 0.0000011*cos2 + 0.0000002*cos3 + ...
        0.0000033*cos4) + 0.000000208);
    Gnm(10,11,:) = -(0.96*(0.0000025*cos1 - 0.000003*cos2 + 0.0000005*cos3 + 0.0000009*cos4) - ...
        0.00000115);

end