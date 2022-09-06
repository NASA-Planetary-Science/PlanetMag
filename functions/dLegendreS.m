% Returns the value of the derivative of the Schmidt-seminormalized Legendre function of degree n and order m given the angle theta
function dPnm = dLegendreS(n, m, theta, UNNORM)
    if ~exist('UNNORM', 'var'); UNNORM = 0; end
    
    if (n == 1)
        if     (m == 0); dPnm = -sin(theta);
        elseif (m == 1); dPnm = cos(theta);
        else;            dPnm = 0;
        end
    elseif (n == 2)
        if     (m == 0); dPnm = -3 * cos(theta) .* sin(theta);
        elseif (m == 1); dPnm = sqrt(3) * cos(2*theta);
        elseif (m == 2); dPnm = sqrt(3) * cos(theta) .* sin(theta);
        else;            dPnm = 0;
        end
    elseif (n == 3)
        if     (m == 0); dPnm = -(3/8) * (sin(theta) + 5*sin(3*theta));
        elseif (m == 1); dPnm = sqrt(3/2)/8 * (cos(theta) + 15*cos(3*theta));    
        elseif (m == 2); dPnm = -sqrt(15)/8 * (sin(theta) - 3*sin(3*theta)); 
        elseif (m == 3); dPnm = sqrt(5/2)*3/2 * cos(theta) .* sin(theta).^2;        
        else;            dPnm = 0;
        end
    elseif (n == 4)
        if     (m == 0); dPnm = -(5/2) * cos(theta) .* (7*cos(theta).^2 - 3) .* sin(theta);
        elseif (m == 1); dPnm = sqrt(10)/4 * (28*cos(theta).^4 - 27*cos(theta).^2 + 3);
        elseif (m == 2); dPnm = sqrt(5) * cos(theta) .* (7*cos(theta).^2 - 4) .* sin(theta);
        elseif (m == 3); dPnm = sqrt(70)/4 * (2*cos(theta).^2 - 1) .* sin(theta).^2;  
        elseif (m == 4); dPnm = sqrt(35)/2 * cos(theta) .* sin(theta).^3;
        else;            dPnm = 0;
        end
    elseif (n == 5)
        if     (m == 0); dPnm = -(15/8) * (21*cos(theta).^4 - 14*cos(theta).^2 + 1) .* sin(theta);
        elseif (m == 1); dPnm = sqrt(15)/8 * cos(theta) .* (105*cos(theta).^4 - 126*cos(theta).^2 + 29);
        elseif (m == 2); dPnm = sqrt(105)/4 * (15*cos(theta).^4 - 12*cos(theta).^2 + 1) .* sin(theta);
        elseif (m == 3); dPnm = sqrt(35/2)*3/8 * cos(theta) .* (15*cos(theta).^2 - 7) .* sin(theta).^2;
        elseif (m == 4); dPnm = -sqrt(35)*3/8 * (5*sin(theta).^2 - 4) .* sin(theta).^3;
        elseif (m == 5); dPnm = sqrt(7/2)*15/8 * cos(theta) .* sin(theta).^4;
        else;            dPnm = 0;
        end
    elseif (n == 6)
        if     (m == 0); dPnm = -(21/8) * cos(theta) .* (33*cos(theta).^4 - 30*cos(theta).^2 + 5) .* sin(theta);
        elseif (m == 1); dPnm = sqrt(21)/8 * (198*cos(theta).^6 - 285*cos(theta).^4 + 100*cos(theta).^2 - 5);
        elseif (m == 2); dPnm = sqrt(105/2)/16 * cos(theta) .* (198*cos(theta).^4 - 204*cos(theta).^2 + 38) .* sin(theta);
        elseif (m == 3); dPnm = sqrt(105/2)/8 * (66*sin(theta).^4 - 87*sin(theta).^2 + 24) .* sin(theta).^2;
        elseif (m == 4); dPnm = sqrt(7)*3/8 * cos(theta) .* (33*cos(theta).^2 - 13) .* sin(theta).^3;
        elseif (m == 5); dPnm = -sqrt(77/2)*3/8 * (6*sin(theta).^2 - 5) .* sin(theta).^4;
        elseif (m == 6); dPnm = sqrt(231/2)*3/8 * cos(theta) .* sin(theta).^5;
        else;            dPnm = 0;
        end
    elseif (n == 7)
        if     (m == 0); dPnm = -(7/16) * (429*cos(theta).^6 - 495*cos(theta).^4 + 135*cos(theta).^2 - 5) .* sin(theta);
        elseif (m == 1); dPnm = sqrt(7)/32 * cos(theta) .* (3003*cos(theta).^6 - 5049*cos(theta).^4 + 2385*cos(theta).^2 - 275);
        elseif (m == 2); dPnm = sqrt(21/2)/16 * (1001*cos(theta).^6 - 1265*cos(theta).^4 + 375*cos(theta).^2 - 15) .* sin(theta);
        elseif (m == 3); dPnm = sqrt(21)/32 * cos(theta) .* (1001*cos(theta).^4 - 902*cos(theta).^2 + 141) .* sin(theta).^2;
        elseif (m == 4); dPnm = sqrt(231)/16 * (91*cos(theta).^4 - 54*cos(theta).^2 + 3) .* sin(theta).^3;
        elseif (m == 5); dPnm = sqrt(231)/32 * cos(theta) .* (91*cos(theta).^2 - 31) .* sin(theta).^4;
        elseif (m == 6); dPnm = sqrt(3003/2)/16 * (7*cos(theta).^2 - 1) .* sin(theta).^5;
        elseif (m == 7); dPnm = sqrt(429)*7/32 * cos(theta) .* sin(theta).^6;
        else;            dPnm = 0;
        end
    elseif (n == 8)
        if     (m == 0); dPnm = -(9/16) * cos(theta) .* (715*cos(theta).^6 - 1001*cos(theta).^4 + 385*cos(theta).^2 - 35) .* sin(theta);
        elseif (m == 1); dPnm = (3/32) * (5720*cos(theta).^8 - 11011*cos(theta).^6 + 6545*cos(theta).^4 - 1225*cos(theta).^2 + 35);
        elseif (m == 2); dPnm = sqrt(35/2)*3/8 * cos(theta) .* (286*cos(theta).^6 - 429*cos(theta).^4 + 176*cos(theta).^2 - 17) .* sin(theta);
        elseif (m == 3); dPnm = sqrt(1155)*3/32 * (104*cos(theta).^6 - 117*cos(theta).^4 + 30*cos(theta).^2 - 1) .* sin(theta).^2;
        elseif (m == 4); dPnm = sqrt(77)*3/8 * cos(theta) .* (65*cos(theta).^4 - 52*cos(theta).^2 + 7) .* sin(theta).^3;
        elseif (m == 5); dPnm = sqrt(1001)*3/32 * (40*cos(theta).^4 - 21*cos(theta).^2 + 1) .* sin(theta).^4;
        elseif (m == 6); dPnm = sqrt(429/2)*3/8 * cos(theta) .* (10*cos(theta).^2 - 3) .* sin(theta).^5;
        elseif (m == 7); dPnm = sqrt(715)*3/32 * (8*cos(theta).^2 - 1) .* sin(theta).^6;
        elseif (m == 8); dPnm = sqrt(715)*3/16 * cos(theta) .* sin(theta).^7;
        else;            dPnm = 0;
        end
    elseif (n == 9)
        if     (m == 0); dPnm = -(45/128) * (2431*cos(theta).^8 - 4004*cos(theta).^6 + 2002*cos(theta).^4 - 308*cos(theta).^2 + 7) .* sin(theta);
        elseif (m == 1); dPnm = sqrt(5)*3/128 * cos(theta) .* (21879*cos(theta).^8 - 47476*cos(theta).^6 + 34034*cos(theta).^4 - 8932*cos(theta).^2 + 623);
        elseif (m == 2); dPnm = sqrt(55/2)*3/32 * (1989*cos(theta).^8 - 3458*cos(theta).^6 + 1820*cos(theta).^4 - 294*cos(theta).^2 + 7) .* sin(theta);
        elseif (m == 3); dPnm = sqrt(1155/2)*3/64 * cos(theta) .* (663*cos(theta).^6 - 897*cos(theta).^4 + 325*cos(theta).^2 - 27) .* sin(theta).^2;
        elseif (m == 4); dPnm = sqrt(5005)*3/64 * (153*cos(theta).^6 - 155*cos(theta).^4 + 35*cos(theta).^2 - 1) .* sin(theta).^3;
        elseif (m == 5); dPnm = sqrt(143/2)*15/64 * cos(theta) .* (153*cos(theta).^4 - 110*cos(theta).^2 + 13) .* sin(theta).^4;
        elseif (m == 6); dPnm = sqrt(2145/2)*3/32 * (51*cos(theta).^4 - 24*cos(theta).^2 + 1) .* sin(theta).^5;
        elseif (m == 7); dPnm = sqrt(715/2)*3/128 * cos(theta) .* (153*cos(theta).^2 - 41) .* sin(theta).^6;
        elseif (m == 8); dPnm = sqrt(12155)*3/128 * (9*cos(theta).^2 - 1) .* sin(theta).^7;
        elseif (m == 9); dPnm = sqrt(12155/2)*9/128 * cos(theta) .* sin(theta).^8;
        else;            dPnm = 0;
        end
    elseif (n == 10)
        if     (m == 0);  dPnm = -(55/128) * cos(theta) .* (4199*cos(theta).^8 - 7956*cos(theta).^6 + 4914*cos(theta).^4 - 1092*cos(theta).^2 + 63) .* sin(theta);
        elseif (m == 1);  dPnm = sqrt(55)/128 * (41990*cos(theta).^10 - 101439*cos(theta).^8 + 85176*cos(theta).^6 - 28938*cos(theta).^4 + 3402*cos(theta).^2 - 63);
        elseif (m == 2);  dPnm = sqrt(165)/128 * cos(theta) .* (20995*cos(theta).^8 - 41548*cos(theta).^6 + 26754*cos(theta).^4 - 6188*cos(theta).^2 + 371) .* sin(theta);
        elseif (m == 3);  dPnm = sqrt(2145/2)/64 * (3230*cos(theta).^8 - 5117*cos(theta).^6 + 2415*cos(theta).^4 - 343*cos(theta).^2 + 7) .* sin(theta).^2;
        elseif (m == 4);  dPnm = sqrt(2145)/64 * cos(theta) .* (1615*cos(theta).^6 - 1989*cos(theta).^4 + 645*cos(theta).^2 - 47) .* sin(theta).^3;
        elseif (m == 5);  dPnm = sqrt(429/2)*5/64 * (646*cos(theta).^6 - 595*cos(theta).^4 + 120*cos(theta).^2 - 3) .* sin(theta).^4;
        elseif (m == 6);  dPnm = sqrt(2145/2)/128 * cos(theta) .* (1615*cos(theta).^4 - 1054*cos(theta).^2 + 111) .* sin(theta).^5;
        elseif (m == 7);  dPnm = sqrt(36465/2)/128 * (190*cos(theta).^4 - 81*cos(theta).^2 + 3) .* sin(theta).^6;
        elseif (m == 8);  dPnm = sqrt(12155)/128 * cos(theta) .* (95*cos(theta).^2 - 23) .* sin(theta).^7;
        elseif (m == 9);  dPnm = sqrt(230945/2)/128 * (10*cos(theta).^2 - 1) .* sin(theta).^8;
        elseif (m == 10); dPnm = sqrt(46189/2)*5/128 * cos(theta) .* sin(theta).^9;
        else;             dPnm = 0;
        end
    else; dPnm = 0;
    end
    
    if UNNORM && m ~= 0
        SchmidtNormFac = sqrt(2*factorial(n-m)/factorial(n+m));
        dPnm = dPnm / SchmidtNormFac;
    end
end