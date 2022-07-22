% Returns the value of Schmidt-seminormalized Legendre function of degree n and order m given the angle theta
function Pnm = LegendreS(n, m, theta, UNNORM)  % n ( or l) = degree m = order
    costh = cos(theta);
    sinth = sin(theta);
    if ~exist('UNNORM', 'var'); UNNORM = 0; end

    if (n == 1)
        if     (m == 0); Pnm = costh;
        elseif (m == 1); Pnm = sinth;
        else;            Pnm = 0;
        end
    elseif (n == 2)
        if     (m == 0); Pnm = (1/2) * (3*costh.^2 - 1);
        elseif (m == 1); Pnm = sqrt(3) * costh .* sinth;
        elseif (m == 2); Pnm = sqrt(3)/2 * sinth.^2;
        else;            Pnm = 0;
        end
    elseif (n == 3)
        if     (m == 0); Pnm = (1/2) * (5*costh.^3 - 3*costh);
        elseif (m == 1); Pnm = sqrt(6)/4 * (5*costh.^2 - 1) .* sinth;
        elseif (m == 2); Pnm = sqrt(15)/2 * costh .* sinth.^2;
        elseif (m == 3); Pnm = sqrt(10)/4 * sinth.^3;
        else;            Pnm = 0;
        end
    elseif (n == 4)
        if     (m == 0); Pnm = (1/8) * (35*costh.^4 - 30*costh.^2 + 3);
        elseif (m == 1); Pnm = sqrt(10)/4 * (7*costh.^3 - 3*costh) .* sinth;
        elseif (m == 2); Pnm = sqrt(5)/4 * (7*costh.^2 - 1) .* sinth.^2;
        elseif (m == 3); Pnm = sqrt(70)/4 * costh .* sinth.^3;
        elseif (m == 4); Pnm = sqrt(35)/8 * sinth.^4;
        else;            Pnm = 0;
        end
    elseif (n == 5)
        if     (m == 0); Pnm = (1/8) * costh .* (63*costh.^4 - 70*costh.^2 + 15);
        elseif (m == 1); Pnm = sqrt(15)/8 * (21*costh.^4 - 14*costh.^2 + 1) .* sinth;
        elseif (m == 2); Pnm = sqrt(105)/4 * costh .* (3*costh.^2 - 1) .* sinth.^2;
        elseif (m == 3); Pnm = sqrt(35/2)/8 * (9*costh.^2 - 1) .* sinth.^3;
        elseif (m == 4); Pnm = sqrt(35)*3/8 * costh .* sinth.^4;
        elseif (m == 5); Pnm = sqrt(7/2)*3/8 * sinth.^5;
        else;            Pnm = 0;
        end
    elseif (n == 6)
        if     (m == 0); Pnm = (1/16) * (231*costh.^6 - 315*costh.^4 + 105*costh.^2 - 5);
        elseif (m == 1); Pnm = sqrt(21)/8 * costh .* (33*costh.^4 - 30*costh.^2 + 5) .* sinth;
        elseif (m == 2); Pnm = sqrt(105/2)/16 * (33*costh.^4 - 18*costh.^2 + 1) .* sinth.^2;
        elseif (m == 3); Pnm = sqrt(105/2)/8 * costh .* (11*costh.^2 - 3) .* sinth.^3;
        elseif (m == 4); Pnm = sqrt(7)*3/16 * (11*costh.^2 - 1) .* sinth.^4;
        elseif (m == 5); Pnm = sqrt(77/2)*3/8 * costh .* sinth.^5;
        elseif (m == 6); Pnm = sqrt(231/2)/16 * sinth.^6;
        else;            Pnm = 0;
        end
    elseif (n == 7)
        if     (m == 0); Pnm = (1/16) * (429*costh.^7 - 693*costh.^5 + 315*costh.^3 - 35*costh);
        elseif (m == 1); Pnm = sqrt(7)/32 * (429*costh.^6 - 495*costh.^4 + 135*costh.^2 - 5) .* sinth;
        elseif (m == 2); Pnm = sqrt(21/2)/16 * costh .* (143*costh.^4 - 110*costh.^2 + 15) .* sinth.^2;
        elseif (m == 3); Pnm = sqrt(21)/32 * (143*costh.^4 - 66*costh.^2 + 3) .* sinth.^3;
        elseif (m == 4); Pnm = sqrt(231)/16 * costh .* (13*costh.^2 - 3) .* sinth.^4;
        elseif (m == 5); Pnm = sqrt(231)/32 * (13*costh.^2 - 1) .* sinth.^5;
        elseif (m == 6); Pnm = sqrt(3003/2)/16 * costh .* sinth.^6;
        elseif (m == 7); Pnm = sqrt(429)/32 * sinth.^7;
        else;            Pnm = 0;
        end
    elseif (n == 8)
        if     (m == 0); Pnm = (1/128) * (6435*costh.^8 - 12012*costh.^6 + 6930*costh.^4 - 1260*costh.^2 + 35);
        elseif (m == 1); Pnm = 3/32 * costh .* (715*costh.^6 - 1001*costh.^4 + 385*costh.^2 - 35) .* sinth;
        elseif (m == 2); Pnm = sqrt(35/2)*3/32 * (143*costh.^6 - 143*costh.^4 + 33*costh.^2 - 1) .* sinth.^2;
        elseif (m == 3); Pnm = sqrt(1155)/32 * costh .* (39*costh.^4 - 26*costh.^2 + 3) .* sinth.^3;
        elseif (m == 4); Pnm = sqrt(77)*3/64 * (65*costh.^4 - 26*costh.^2 + 1) .* sinth.^4;
        elseif (m == 5); Pnm = sqrt(1001)*3/32 * costh .* (5*costh.^2 - 1) .* sinth.^5;
        elseif (m == 6); Pnm = sqrt(429/2)/32 * (15*costh.^2 - 1) .* sinth.^6;
        elseif (m == 7); Pnm = sqrt(715)*3/32 * costh .* sinth.^7;
        elseif (m == 8); Pnm = sqrt(715)*3/128 * sinth.^8;
        else;            Pnm = 0;
        end
    elseif (n == 9)
        if     (m == 0); Pnm = (1/128) * costh .* (12155*costh.^8 - 25740*costh.^6 + 18018*costh.^4 - 4620*costh.^2 + 315);
        elseif (m == 1); Pnm = sqrt(5)*3/128 * (2431*costh.^8 - 4004*costh.^6 + 2002*costh.^4 - 308*costh.^2 + 7) .* sinth;
        elseif (m == 2); Pnm = sqrt(55/2)*3/32 * costh .* (221*costh.^6 - 273*costh.^4 + 91*costh.^2 - 7) .* sinth.^2;
        elseif (m == 3); Pnm = sqrt(1155/2)/64 * (221*costh.^6 - 195*costh.^4 + 39*costh.^2 - 1) .* sinth.^3;
        elseif (m == 4); Pnm = sqrt(5005)*3/64 * costh .* (17*costh.^4 - 10*costh.^2 + 1) .* sinth.^4;
        elseif (m == 5); Pnm = sqrt(143/2)*3/64 * (85*costh.^4 - 30*costh.^2 + 1) .* sinth.^5;
        elseif (m == 6); Pnm = sqrt(2145/2)/32 * costh .* (17*costh.^2 - 3) .* sinth.^6;
        elseif (m == 7); Pnm = sqrt(715/2)*3/128 * (17*costh.^2 - 1) .* sinth.^7;
        elseif (m == 8); Pnm = sqrt(12155)*3/128 * costh .* sinth.^8;
        elseif (m == 9); Pnm = sqrt(12155/2)/128 * sinth.^9;
        else;            Pnm = 0;
        end
    elseif (n == 10)
        if     (m == 0);  Pnm = (1/256) * (46189*costh.^10 - 109395*costh.^8 + 90090*costh.^6 - 30030*costh.^4 + 3465*costh.^2 - 63);
        elseif (m == 1);  Pnm = sqrt(55)/128 * costh .* (4199*costh.^8 - 7956*costh.^6 + 4914*costh.^4 - 1092*costh.^2 + 63) .* sinth;
        elseif (m == 2);  Pnm = sqrt(165)/256 * (4199*costh.^8 - 6188*costh.^6 + 2730*costh.^4 - 364*costh.^2 + 7) .* sinth.^2;
        elseif (m == 3);  Pnm = sqrt(2145/2)/64 * costh .* (323*costh.^6 - 357*costh.^4 + 105*costh.^2 - 7) .* sinth.^3;
        elseif (m == 4);  Pnm = sqrt(2145)/128 * (323*costh.^6 - 255*costh.^4 + 45*costh.^2 - 1) .* sinth.^4;
        elseif (m == 5);  Pnm = sqrt(429/2)/64 * costh .* (323*costh.^4 - 170*costh.^2 + 15) .* sinth.^5;
        elseif (m == 6);  Pnm = sqrt(2145/2)/256 * (323*costh.^6 - 102*costh.^2 + 3) .* sinth.^6;
        elseif (m == 7);  Pnm = sqrt(36465/2)/128 * costh .* (19*costh.^2 - 3) .* sinth.^7;
        elseif (m == 8);  Pnm = sqrt(12155)/256 * (19*costh.^2 - 1) .* sinth.^8;
        elseif (m == 9);  Pnm = sqrt(230945/2)/128 * costh .* sinth.^9;
        elseif (m == 10); Pnm = sqrt(46189/2)/256 * sinth.^10;
        else;             Pnm = 0;
        end
    else; Pnm = 0;
    end
    
    if UNNORM && m ~= 0
        SchmidtNormFac = sqrt(2*factorial(n+m)/factorial(n-m));
        Pnm = Pnm / SchmidtNormFac;
    end
end