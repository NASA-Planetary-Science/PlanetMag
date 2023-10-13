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
        if     (m == 0); Pnm = (1/64) * (9 + 20*cos(2*theta) + 35*cos(4*theta));
        elseif (m == 1); Pnm = sqrt(5/2)/8 * (9*cos(theta) + 7*cos(3*theta)) .* sin(theta);
        elseif (m == 2); Pnm = sqrt(5)/8 * (5 + 7*cos(2*theta)) .* sin(theta).^2;
        elseif (m == 3); Pnm = sqrt(35/2)/2 * cos(theta) .* sin(theta).^3;
        elseif (m == 4); Pnm = sqrt(35)/8 * sin(theta).^4;
        else;            Pnm = 0;
        end
    elseif (n == 5)
        if     (m == 0); Pnm = (1/128) * (30*cos(theta) + 35*cos(3*theta) + 63*cos(5*theta));
        elseif (m == 1); Pnm = sqrt(15)/64 * (15 + 28*cos(2*theta) + 21*cos(4*theta)) .* sin(theta);
        elseif (m == 2); Pnm = sqrt(105)/16 * (5*cos(theta) + 3*cos(3*theta)) .* sin(theta).^2;
        elseif (m == 3); Pnm = sqrt(35/2)/16 * (7 + 9*cos(2*theta)) .* sin(theta).^3;
        elseif (m == 4); Pnm = sqrt(35)*3/8 * cos(theta) .* sin(theta).^4;
        elseif (m == 5); Pnm = sqrt(7/2)*3/8 * sin(theta).^5;
        else;            Pnm = 0;
        end
    elseif (n == 6)
        if     (m == 0); Pnm = (1/512) * (50 + 105*cos(2*theta) + 126*cos(4*theta) + 231*cos(6*theta));
        elseif (m == 1); Pnm = sqrt(21)/128 * (50*cos(theta) + 45*cos(3*theta) + 33*cos(5*theta)) .* sin(theta);
        elseif (m == 2); Pnm = sqrt(105/2)/128 * (35 + 60*cos(2*theta) + 33*cos(4*theta)) .* sin(theta).^2;
        elseif (m == 3); Pnm = sqrt(105/2)/32 * (21*cos(theta) + 11*cos(3*theta)) .* sin(theta).^3;
        elseif (m == 4); Pnm = sqrt(7)*3/32 * (9 + 11*cos(2*theta)) .* sin(theta).^4;
        elseif (m == 5); Pnm = sqrt(77/2)*3/8 * cos(theta) .* sin(theta).^5;
        elseif (m == 6); Pnm = sqrt(231/2)/16 * sin(theta).^6;
        else;            Pnm = 0;
        end
    elseif (n == 7)
        if     (m == 0); Pnm = (1/1024) * (175*cos(theta) + 189*cos(3*theta) + 231*cos(5*theta) + 429*cos(7*theta));
        elseif (m == 1); Pnm = sqrt(7)/1024 * (350 + 675*cos(2*theta) + 594*cos(4*theta) + 429*cos(6*theta)) .* sin(theta);
        elseif (m == 2); Pnm = sqrt(21/2)/256 * (350*cos(theta) + 275*cos(3*theta) + 143*cos(5*theta)) .* sin(theta).^2;
        elseif (m == 3); Pnm = sqrt(21)/256 * (189 + 308*cos(2*theta) + 143*cos(4*theta)) .* sin(theta).^3;
        elseif (m == 4); Pnm = sqrt(231)/64 * (27*cos(theta) + 13*cos(3*theta)) .* sin(theta).^4;
        elseif (m == 5); Pnm = sqrt(231)/64 * (11 + 13*cos(2*theta)) .* sin(theta).^5;
        elseif (m == 6); Pnm = sqrt(3003/2)/16 * cos(theta) .* sin(theta).^6;
        elseif (m == 7); Pnm = sqrt(429)/32 * sin(theta).^7;
        else;            Pnm = 0;
        end
    elseif (n == 8)
        if     (m == 0); Pnm = (1/16384) * (1225 + 2520*cos(2*theta) + 2772*cos(4*theta) + 3432*cos(6*theta) + 6435*cos(8*theta));
        elseif (m == 1); Pnm = (3/2048) * (1225*cos(theta) + 11*(105*cos(3*theta) + 91*cos(5*theta) + 65*cos(7*theta))) .* sin(theta);
        elseif (m == 2); Pnm = sqrt(35/2)*3/1024 * (210 + 385*cos(2*theta) + 286*cos(4*theta) + 143*cos(6*theta)) .* sin(theta).^2;
        elseif (m == 3); Pnm = sqrt(1155)/512 * (126*cos(theta) + 91*cos(3*theta) + 39*cos(5*theta)) .* sin(theta).^3;
        elseif (m == 4); Pnm = sqrt(77)*3/512 * (99 + 156*cos(2*theta) + 65*cos(4*theta)) .* sin(theta).^4;
        elseif (m == 5); Pnm = sqrt(1001)*3/128 * (11*cos(theta) + 5*cos(3*theta)) .* sin(theta).^5;
        elseif (m == 6); Pnm = sqrt(429/2)/64 * (13 + 15*cos(2*theta)) .* sin(theta).^6;
        elseif (m == 7); Pnm = sqrt(715)*3/32 * cos(theta) .* sin(theta).^7;
        elseif (m == 8); Pnm = sqrt(715)*3/128 * sin(theta).^8;
        else;            Pnm = 0;
        end
    elseif (n == 9)
        if     (m == 0); Pnm = (1/32768) * (4410*cos(theta) + 4620*cos(3*theta) + 143*(36*cos(5*theta) + 45*cos(7*theta) + 85*cos(9*theta)));
        elseif (m == 1); Pnm = sqrt(5)*3/16384 * (2205 + 4312*cos(2*theta) + 4004*cos(4*theta) + 3432*cos(6*theta) + 2431*cos(8*theta)) .* sin(theta);
        elseif (m == 2); Pnm = sqrt(55/2)*3/2048 * (735*cos(theta) + 637*cos(3*theta) + 455*cos(5*theta) + 221*cos(7*theta)) .* sin(theta).^2;
        elseif (m == 3); Pnm = sqrt(1155/2)/2048 * (462 + 819*cos(2*theta) + 546*cos(4*theta) + 221*cos(6*theta)) .* sin(theta).^3;
        elseif (m == 4); Pnm = sqrt(5005)*3/1024 * (66*cos(theta) + 45*cos(3*theta) + 17*cos(5*theta)) .* sin(theta).^4;
        elseif (m == 5); Pnm = sqrt(143/2)*3/512 * (143 + 220*cos(2*theta) + 85*cos(4*theta)) .* sin(theta).^5;
        elseif (m == 6); Pnm = sqrt(2145/2)/128 * (39*cos(theta) + 17*cos(3*theta)) .* sin(theta).^6;
        elseif (m == 7); Pnm = sqrt(715/2)*3/256 * (15 + 17*cos(2*theta)) .* sin(theta).^7;
        elseif (m == 8); Pnm = sqrt(12155)*3/128 * cos(theta) .* sin(theta).^8;
        elseif (m == 9); Pnm = sqrt(12155/2)/128 * sin(theta).^9;
        else;            Pnm = 0;
        end
    elseif (n == 10)
        if     (m == 0);  Pnm = (1/131072) * (7938 + 16170*cos(2*theta) + 17160*cos(4*theta) + 19305*cos(6*theta) + 24310*cos(8*theta) + 46189*cos(10*theta));
        elseif (m == 1);  Pnm = sqrt(55)/32768 * (7938*cos(theta) + 13*(588*cos(3*theta) + 540*cos(5*theta) + 459*cos(7*theta) + 323*cos(9*theta))) .* sin(theta);
        elseif (m == 2);  Pnm = sqrt(165)/32768 * (8085 + 15288*cos(2*theta) + 12740*cos(4*theta) + 8840*cos(6*theta) + 4199*cos(8*theta)) .* sin(theta).^2;
        elseif (m == 3);  Pnm = sqrt(2145/2)/4096 * (1617*cos(theta) + 1323*cos(3*theta) + 833*cos(5*theta) + 323*cos(7*theta)) .* sin(theta).^3;
        elseif (m == 4);  Pnm = sqrt(2145)/4096 * (858 + 1485*cos(2*theta) + 918*cos(4*theta) + 323*cos(6*theta)) .* sin(theta).^4;
        elseif (m == 5);  Pnm = sqrt(429/2)/1024 * (1430*cos(theta) + 935*cos(3*theta) + 323*cos(5*theta)) .* sin(theta).^5;
        elseif (m == 6);  Pnm = sqrt(2145/2)/2048 * (585 + 884*cos(2*theta) + 323*cos(4*theta)) .* sin(theta).^6;
        elseif (m == 7);  Pnm = sqrt(36465/2)/512 * (45*cos(theta) + 19*cos(3*theta)) .* sin(theta).^7;
        elseif (m == 8);  Pnm = sqrt(12155)/512 * (17 + 19*cos(2*theta)) .* sin(theta).^8;
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