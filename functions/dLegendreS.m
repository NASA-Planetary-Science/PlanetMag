% Returns the value of the derivative of the Schmidt-seminormalized Legendre function of degree n and order m given the angle theta
function dPnm = dLegendreS(n, m, theta, UNNORM)
    costh = cos(theta);
    sinth = sin(theta);
    if ~exist('UNNORM', 'var'); UNNORM = 0; end
    
    if (n == 1)
        if     (m == 0); dPnm = -sinth;
        elseif (m == 1); dPnm = costh;
        else;            dPnm = 0;
        end
    elseif (n == 2)
        if     (m == 0); dPnm = -3 * costh .* sinth;
        elseif (m == 1); dPnm = sqrt(3) * (2*costh.^2 - 1);
        elseif (m == 2); dPnm = sqrt(3) * costh .* sinth;
        else;            dPnm = 0;
        end
    elseif (n == 3)
        if     (m == 0); dPnm = -(3/2) * (10*costh/3 - 1) .* sinth;
        elseif (m == 1); dPnm = -sqrt(6)/4 * (15*sinth.^2 - 4) .* costh;    
        elseif (m == 2); dPnm = sqrt(15)/2 * (3*costh.^2 - 1) .* sinth; 
        elseif (m == 3); dPnm = sqrt(10)*3/4 * costh .* sinth.^2;        
        else;            dPnm = 0;
        end
    elseif (n == 4)
        if     (m == 0); dPnm = -(5/2) * costh .* (7*costh.^2 - 3) .* sinth;
        elseif (m == 1); dPnm = sqrt(10)/4 * (28*costh.^4 - 27*costh.^2 + 3);
        elseif (m == 2); dPnm = sqrt(5) * costh .* (7*costh.^2 - 4) .* sinth;
        elseif (m == 3); dPnm = sqrt(70)/4 * (2*costh.^2 - 1) .* sinth.^2;  
        elseif (m == 4); dPnm = sqrt(35)/2 * costh .* sinth.^3;
        else;            dPnm = 0;
        end
    elseif (n == 5)
        if     (m == 0); dPnm = -(15/8) * (21*costh.^4 - 14*costh.^2 + 1) .* sinth;
        elseif (m == 1); dPnm = sqrt(15)/8 * costh .* (105*costh.^4 - 126*costh.^2 + 29);
        elseif (m == 2); dPnm = sqrt(105)/4 * (15*costh.^4 - 12*costh.^2 + 1) .* sinth;
        elseif (m == 3); dPnm = sqrt(35/2)*3/8 * costh .* (15*costh.^2 - 7) .* sinth.^2;
        elseif (m == 4); dPnm = -sqrt(35)*3/8 * (5*sinth.^2 - 4) .* sinth.^3;
        elseif (m == 5); dPnm = sqrt(7/2)*15/8 * costh .* sinth.^4;
        else;            dPnm = 0;
        end
    elseif (n == 6)
        if     (m == 0); dPnm = -(21/8) * costh .* (33*costh.^4 - 30*costh.^2 + 5) .* sinth;
        elseif (m == 1); dPnm = sqrt(21)/8 * (198*costh.^6 - 285*costh.^4 + 100*costh.^2 - 5);
        elseif (m == 2); dPnm = sqrt(105/2)/16 * costh .* (198*costh.^4 - 204*costh.^2 + 38) .* sinth;
        elseif (m == 3); dPnm = sqrt(105/2)/8 * (66*sinth.^4 - 87*sinth.^2 + 24) .* sinth.^2;
        elseif (m == 4); dPnm = sqrt(7)*3/8 * costh .* (33*costh.^2 - 13) .* sinth.^3;
        elseif (m == 5); dPnm = -sqrt(77/2)*3/8 * (6*sinth.^2 - 5) .* sinth.^4;
        elseif (m == 6); dPnm = sqrt(231/2)*3/8 * costh .* sinth.^5;
        else;            dPnm = 0;
        end
    elseif (n == 7)
        if     (m == 0); dPnm = -(7/16) * (429*costh.^6 - 495*costh.^4 + 135*costh.^2 - 5) .* sinth;
        elseif (m == 1); dPnm = sqrt(7)/32 * costh .* (3003*costh.^6 - 5049*costh.^4 + 2385*costh.^2 - 275);
        elseif (m == 2); dPnm = sqrt(21/2)/16 * (1001*costh.^6 - 1265*costh.^4 + 375*costh.^2 - 15) .* sinth;
        elseif (m == 3); dPnm = sqrt(21)/32 * costh .* (1001*costh.^4 - 902*costh.^2 + 141) .* sinth.^2;
        elseif (m == 4); dPnm = sqrt(231)/16 * (91*costh.^4 - 54*costh.^2 + 3) .* sinth.^3;
        elseif (m == 5); dPnm = sqrt(231)/32 * costh .* (91*costh.^2 - 31) .* sinth.^4;
        elseif (m == 6); dPnm = sqrt(3003/2)/16 * (7*costh.^2 - 1) .* sinth.^5;
        elseif (m == 7); dPnm = sqrt(429)*7/32 * costh .* sinth.^6;
        else;            dPnm = 0;
        end
    elseif (n == 8)
        if     (m == 0); dPnm = -(9/16) * costh .* (715*costh.^6 - 1001*costh.^4 + 385*costh.^2 - 35) .* sinth;
        elseif (m == 1); dPnm = (3/32) * (5720*costh.^8 - 11011*costh.^6 + 6545*costh.^4 - 1225*costh.^2 + 35);
        elseif (m == 2); dPnm = sqrt(35/2)*3/8 * costh .* (286*costh.^6 - 429*costh.^4 + 176*costh.^2 - 17) .* sinth;
        elseif (m == 3); dPnm = sqrt(1155)*3/32 * (104*costh.^6 - 117*costh.^4 + 30*costh.^2 - 1) .* sinth.^2;
        elseif (m == 4); dPnm = sqrt(77)*3/8 * costh .* (65*costh.^4 - 52*costh.^2 + 7) .* sinth.^3;
        elseif (m == 5); dPnm = sqrt(1001)*3/32 * (40*costh.^4 - 21*costh.^2 + 1) .* sinth.^4;
        elseif (m == 6); dPnm = sqrt(429/2)*3/8 * costh .* (10*costh.^2 - 3) .* sinth.^5;
        elseif (m == 7); dPnm = sqrt(715)*3/32 * (8*costh.^2 - 1) .* sinth.^6;
        elseif (m == 8); dPnm = sqrt(715)*3/16 * costh .* sinth.^7;
        else;            dPnm = 0;
        end
    elseif (n == 9)
        if     (m == 0); dPnm = -(45/128) * (2431*costh.^8 - 4004*costh.^6 + 2002*costh.^4 - 308*costh.^2 + 7) .* sinth;
        elseif (m == 1); dPnm = sqrt(5)*3/128 * costh .* (21879*costh.^8 - 47476*costh.^6 + 34034*costh.^4 - 8932*costh.^2 + 623);
        elseif (m == 2); dPnm = sqrt(55/2)*3/32 * (1989*costh.^8 - 3458*costh.^6 + 1820*costh.^4 - 294*costh.^2 + 7) .* sinth;
        elseif (m == 3); dPnm = sqrt(1155/2)*3/64 * costh .* (663*costh.^6 - 897*costh.^4 + 325*costh.^2 - 27) .* sinth.^2;
        elseif (m == 4); dPnm = sqrt(5005)*3/64 * (153*costh.^6 - 155*costh.^4 + 35*costh.^2 - 1) .* sinth.^3;
        elseif (m == 5); dPnm = sqrt(143/2)*15/64 * costh .* (153*costh.^4 - 110*costh.^2 + 13) .* sinth.^4;
        elseif (m == 6); dPnm = sqrt(2145/2)*3/32 * (51*costh.^4 - 24*costh.^2 + 1) .* sinth.^5;
        elseif (m == 7); dPnm = sqrt(715/2)*3/128 * costh .* (153*costh.^2 - 41) .* sinth.^6;
        elseif (m == 8); dPnm = sqrt(12155)*3/128 * (9*costh.^2 - 1) .* sinth.^7;
        elseif (m == 9); dPnm = sqrt(12155/2)*9/128 * costh .* sinth.^8;
        else;            dPnm = 0;
        end
    elseif (n == 10)
        if     (m == 0);  dPnm = -(55/128) * costh .* (4199*costh.^8 - 7956*costh.^6 + 4914*costh.^4 - 1092*costh.^2 + 63) .* sinth;
        elseif (m == 1);  dPnm = sqrt(55)/128 * (41990*costh.^10 - 101439*costh.^8 + 85176*costh.^6 - 28938*costh.^4 + 3402*costh.^2 - 63);
        elseif (m == 2);  dPnm = sqrt(165)/128 * costh .* (20995*costh.^8 - 41548*costh.^6 + 26754*costh.^4 - 6188*costh.^2 + 371) .* sinth;
        elseif (m == 3);  dPnm = sqrt(2145/2)/64 * (3230*costh.^8 - 5117*costh.^6 + 2415*costh.^4 - 343*costh.^2 + 7) .* sinth.^2;
        elseif (m == 4);  dPnm = sqrt(2145)/64 * costh .* (1615*costh.^6 - 1989*costh.^4 + 645*costh.^2 - 47) .* sinth.^3;
        elseif (m == 5);  dPnm = sqrt(429/2)*5/64 * (646*costh.^6 - 595*costh.^4 + 120*costh.^2 - 3) .* sinth.^4;
        elseif (m == 6);  dPnm = sqrt(2145/2)/128 * costh .* (1615*costh.^4 - 1054*costh.^2 + 111) .* sinth.^5;
        elseif (m == 7);  dPnm = sqrt(36465/2)/128 * (190*costh.^4 - 81*costh.^2 + 3) .* sinth.^6;
        elseif (m == 8);  dPnm = sqrt(12155)/128 * costh .* (95*costh.^2 - 23) .* sinth.^7;
        elseif (m == 9);  dPnm = sqrt(230945/2)/128 * (10*costh.^2 - 1) .* sinth.^8;
        elseif (m == 10); dPnm = sqrt(46189/2)*5/128 * costh .* sinth.^9;
        else;             dPnm = 0;
        end
    else; dPnm = 0;
    end
    
    if UNNORM && m ~= 0
        SchmidtNormFac = sqrt(2*factorial(n+m)/factorial(n-m));
        dPnm = dPnm / SchmidtNormFac;
    end
end