% Returns the value of Schmidt-seminormalized Legendre function of degree n and order m given the angle theta
function Pnm = LegendreS(n, m, theta, UNNORM)  % n ( or l) = degree m = order
    if ~exist('UNNORM', 'var'); UNNORM = 0; end

    if (n == 1)
        if     (m == 0); Pnm = cos(theta);
        elseif (m == 1); Pnm = sin(theta);
        else;            Pnm = 0;
        end
    elseif (n == 2)
        if     (m == 0); Pnm = (1/4) * (1 + 3*cos(2*theta));
        elseif (m == 1); Pnm = sqrt(3) * cos(theta) .* sin(theta);
        elseif (m == 2); Pnm = sqrt(3)/2 * sin(theta).^2;
        else;            Pnm = 0;
        end
    elseif (n == 3)
        if     (m == 0); Pnm = (1/8) * (3*cos(theta) + 5*cos(3*theta));
        elseif (m == 1); Pnm = sqrt(3/2)/4 * (3 + 5*cos(2*theta)) .* sin(theta);
        elseif (m == 2); Pnm = sqrt(15)/2 * cos(theta) .* sin(theta).^2;
        elseif (m == 3); Pnm = sqrt(5/2)/2 * sin(theta).^3;
        else;            Pnm = 0;
        end
    elseif (n == 4)
        if     (m == 0); Pnm = (1/8) * (35*cos(theta).^4 - 30*cos(theta).^2 + 3);
        elseif (m == 1); Pnm = sqrt(10)/4 * (7*cos(theta).^3 - 3*cos(theta)) .* sin(theta);
        elseif (m == 2); Pnm = sqrt(5)/4 * (7*cos(theta).^2 - 1) .* sin(theta).^2;
        elseif (m == 3); Pnm = sqrt(70)/4 * cos(theta) .* sin(theta).^3;
        elseif (m == 4); Pnm = sqrt(35)/8 * sin(theta).^4;
        else;            Pnm = 0;
        end
    elseif (n == 5)
        if     (m == 0); Pnm = (1/8) * cos(theta) .* (63*cos(theta).^4 - 70*cos(theta).^2 + 15);
        elseif (m == 1); Pnm = sqrt(15)/8 * (21*cos(theta).^4 - 14*cos(theta).^2 + 1) .* sin(theta);
        elseif (m == 2); Pnm = sqrt(105)/4 * cos(theta) .* (3*cos(theta).^2 - 1) .* sin(theta).^2;
        elseif (m == 3); Pnm = sqrt(35/2)/8 * (9*cos(theta).^2 - 1) .* sin(theta).^3;
        elseif (m == 4); Pnm = sqrt(35)*3/8 * cos(theta) .* sin(theta).^4;
        elseif (m == 5); Pnm = sqrt(7/2)*3/8 * sin(theta).^5;
        else;            Pnm = 0;
        end
    elseif (n == 6)
        if     (m == 0); Pnm = (1/16) * (231*cos(theta).^6 - 315*cos(theta).^4 + 105*cos(theta).^2 - 5);
        elseif (m == 1); Pnm = sqrt(21)/8 * cos(theta) .* (33*cos(theta).^4 - 30*cos(theta).^2 + 5) .* sin(theta);
        elseif (m == 2); Pnm = sqrt(105/2)/16 * (33*cos(theta).^4 - 18*cos(theta).^2 + 1) .* sin(theta).^2;
        elseif (m == 3); Pnm = sqrt(105/2)/8 * cos(theta) .* (11*cos(theta).^2 - 3) .* sin(theta).^3;
        elseif (m == 4); Pnm = sqrt(7)*3/16 * (11*cos(theta).^2 - 1) .* sin(theta).^4;
        elseif (m == 5); Pnm = sqrt(77/2)*3/8 * cos(theta) .* sin(theta).^5;
        elseif (m == 6); Pnm = sqrt(231/2)/16 * sin(theta).^6;
        else;            Pnm = 0;
        end
    elseif (n == 7)
        if     (m == 0); Pnm = (1/16) * (429*cos(theta).^7 - 693*cos(theta).^5 + 315*cos(theta).^3 - 35*cos(theta));
        elseif (m == 1); Pnm = sqrt(7)/32 * (429*cos(theta).^6 - 495*cos(theta).^4 + 135*cos(theta).^2 - 5) .* sin(theta);
        elseif (m == 2); Pnm = sqrt(21/2)/16 * cos(theta) .* (143*cos(theta).^4 - 110*cos(theta).^2 + 15) .* sin(theta).^2;
        elseif (m == 3); Pnm = sqrt(21)/32 * (143*cos(theta).^4 - 66*cos(theta).^2 + 3) .* sin(theta).^3;
        elseif (m == 4); Pnm = sqrt(231)/16 * cos(theta) .* (13*cos(theta).^2 - 3) .* sin(theta).^4;
        elseif (m == 5); Pnm = sqrt(231)/32 * (13*cos(theta).^2 - 1) .* sin(theta).^5;
        elseif (m == 6); Pnm = sqrt(3003/2)/16 * cos(theta) .* sin(theta).^6;
        elseif (m == 7); Pnm = sqrt(429)/32 * sin(theta).^7;
        else;            Pnm = 0;
        end
    elseif (n == 8)
        if     (m == 0); Pnm = (1/128) * (6435*cos(theta).^8 - 12012*cos(theta).^6 + 6930*cos(theta).^4 - 1260*cos(theta).^2 + 35);
        elseif (m == 1); Pnm = 3/32 * cos(theta) .* (715*cos(theta).^6 - 1001*cos(theta).^4 + 385*cos(theta).^2 - 35) .* sin(theta);
        elseif (m == 2); Pnm = sqrt(35/2)*3/32 * (143*cos(theta).^6 - 143*cos(theta).^4 + 33*cos(theta).^2 - 1) .* sin(theta).^2;
        elseif (m == 3); Pnm = sqrt(1155)/32 * cos(theta) .* (39*cos(theta).^4 - 26*cos(theta).^2 + 3) .* sin(theta).^3;
        elseif (m == 4); Pnm = sqrt(77)*3/64 * (65*cos(theta).^4 - 26*cos(theta).^2 + 1) .* sin(theta).^4;
        elseif (m == 5); Pnm = sqrt(1001)*3/32 * cos(theta) .* (5*cos(theta).^2 - 1) .* sin(theta).^5;
        elseif (m == 6); Pnm = sqrt(429/2)/32 * (15*cos(theta).^2 - 1) .* sin(theta).^6;
        elseif (m == 7); Pnm = sqrt(715)*3/32 * cos(theta) .* sin(theta).^7;
        elseif (m == 8); Pnm = sqrt(715)*3/128 * sin(theta).^8;
        else;            Pnm = 0;
        end
    elseif (n == 9)
        if     (m == 0); Pnm = (1/128) * cos(theta) .* (12155*cos(theta).^8 - 25740*cos(theta).^6 + 18018*cos(theta).^4 - 4620*cos(theta).^2 + 315);
        elseif (m == 1); Pnm = sqrt(5)*3/128 * (2431*cos(theta).^8 - 4004*cos(theta).^6 + 2002*cos(theta).^4 - 308*cos(theta).^2 + 7) .* sin(theta);
        elseif (m == 2); Pnm = sqrt(55/2)*3/32 * cos(theta) .* (221*cos(theta).^6 - 273*cos(theta).^4 + 91*cos(theta).^2 - 7) .* sin(theta).^2;
        elseif (m == 3); Pnm = sqrt(1155/2)/64 * (221*cos(theta).^6 - 195*cos(theta).^4 + 39*cos(theta).^2 - 1) .* sin(theta).^3;
        elseif (m == 4); Pnm = sqrt(5005)*3/64 * cos(theta) .* (17*cos(theta).^4 - 10*cos(theta).^2 + 1) .* sin(theta).^4;
        elseif (m == 5); Pnm = sqrt(143/2)*3/64 * (85*cos(theta).^4 - 30*cos(theta).^2 + 1) .* sin(theta).^5;
        elseif (m == 6); Pnm = sqrt(2145/2)/32 * cos(theta) .* (17*cos(theta).^2 - 3) .* sin(theta).^6;
        elseif (m == 7); Pnm = sqrt(715/2)*3/128 * (17*cos(theta).^2 - 1) .* sin(theta).^7;
        elseif (m == 8); Pnm = sqrt(12155)*3/128 * cos(theta) .* sin(theta).^8;
        elseif (m == 9); Pnm = sqrt(12155/2)/128 * sin(theta).^9;
        else;            Pnm = 0;
        end
    elseif (n == 10)
        if     (m == 0);  Pnm = (1/256) * (46189*cos(theta).^10 - 109395*cos(theta).^8 + 90090*cos(theta).^6 - 30030*cos(theta).^4 + 3465*cos(theta).^2 - 63);
        elseif (m == 1);  Pnm = sqrt(55)/128 * cos(theta) .* (4199*cos(theta).^8 - 7956*cos(theta).^6 + 4914*cos(theta).^4 - 1092*cos(theta).^2 + 63) .* sin(theta);
        elseif (m == 2);  Pnm = sqrt(165)/256 * (4199*cos(theta).^8 - 6188*cos(theta).^6 + 2730*cos(theta).^4 - 364*cos(theta).^2 + 7) .* sin(theta).^2;
        elseif (m == 3);  Pnm = sqrt(2145/2)/64 * cos(theta) .* (323*cos(theta).^6 - 357*cos(theta).^4 + 105*cos(theta).^2 - 7) .* sin(theta).^3;
        elseif (m == 4);  Pnm = sqrt(2145)/128 * (323*cos(theta).^6 - 255*cos(theta).^4 + 45*cos(theta).^2 - 1) .* sin(theta).^4;
        elseif (m == 5);  Pnm = sqrt(429/2)/64 * cos(theta) .* (323*cos(theta).^4 - 170*cos(theta).^2 + 15) .* sin(theta).^5;
        elseif (m == 6);  Pnm = sqrt(2145/2)/256 * (323*cos(theta).^6 - 102*cos(theta).^2 + 3) .* sin(theta).^6;
        elseif (m == 7);  Pnm = sqrt(36465/2)/128 * cos(theta) .* (19*cos(theta).^2 - 3) .* sin(theta).^7;
        elseif (m == 8);  Pnm = sqrt(12155)/256 * (19*cos(theta).^2 - 1) .* sin(theta).^8;
        elseif (m == 9);  Pnm = sqrt(230945/2)/128 * cos(theta) .* sin(theta).^9;
        elseif (m == 10); Pnm = sqrt(46189/2)/256 * sin(theta).^10;
        else;             Pnm = 0;
        end
    else; Pnm = 0;
    end
    
    if UNNORM && m ~= 0
        SchmidtNormFac = sqrt(2*factorial(n-m)/factorial(n+m));
        Pnm = Pnm / SchmidtNormFac;
    end
end