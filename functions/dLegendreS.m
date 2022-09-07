% Returns the value of the derivative of the Schmidt-seminormalized Legendre function of degree n and order m given the angle theta
function dPnm = dLegendreS(n, m, theta, UNNORM)
    if ~exist('UNNORM', 'var'); UNNORM = 0; end
    
    if (n == 1)
        if     (m == 0); dPnm = -sin(theta);
        elseif (m == 1); dPnm = cos(theta);
        else;            dPnm = 0;
        end
    elseif (n == 2)
        if     (m == 0); dPnm = (-3) * cos(theta) .* sin(theta);
        elseif (m == 1); dPnm = sqrt(3) * cos(2*theta);
        elseif (m == 2); dPnm = sqrt(3) * cos(theta) .* sin(theta);
        else;            dPnm = 0;
        end
    elseif (n == 3)
        if     (m == 0); dPnm = (-3/8) * (sin(theta) + 5*sin(3*theta));
        elseif (m == 1); dPnm = sqrt(3/2)/8 * (cos(theta) + 15*cos(3*theta));    
        elseif (m == 2); dPnm = -sqrt(15)/8 * (sin(theta) - 3*sin(3*theta)); 
        elseif (m == 3); dPnm = sqrt(5/2)*3/2 * cos(theta) .* sin(theta).^2;        
        else;            dPnm = 0;
        end
    elseif (n == 4)
        if     (m == 0); dPnm = (-5/16) * (2*sin(2*theta) + 7*sin(4*theta));
        elseif (m == 1); dPnm = sqrt(5/2)/4 * (cos(2*theta) + 7*cos(4*theta));
        elseif (m == 2); dPnm = sqrt(5)/8 * (-2*sin(2*theta) + 7*sin(4*theta));
        elseif (m == 3); dPnm = sqrt(35/2)/2 * sin(theta) .* sin(3*theta);  
        elseif (m == 4); dPnm = sqrt(35)/2 * cos(theta) .* sin(theta).^3;
        else;            dPnm = 0;
        end
    elseif (n == 5)
        if     (m == 0); dPnm = (-15/128) * (2*sin(theta) + 7*(sin(3*theta) + 3*sin(5*theta)));
        elseif (m == 1); dPnm = sqrt(15)/128 * (2*cos(theta) + 21*(cos(3*theta) + 5*cos(5*theta)));
        elseif (m == 2); dPnm = sqrt(105)/64 * (-2*sin(theta) - 3*sin(3*theta) + 15*sin(5*theta));
        elseif (m == 3); dPnm = sqrt(35/2)*3/64 * sin(theta) .* (2*sin(2*theta) + 15*sin(4*theta));
        elseif (m == 4); dPnm = sqrt(35)*3/16 * (3 + 5*cos(2*theta)) .* sin(theta).^3;
        elseif (m == 5); dPnm = sqrt(7/2)*15/8 * cos(theta) .* sin(theta).^4;
        else;            dPnm = 0;
        end
    elseif (n == 6)
        if     (m == 0); dPnm = (-21/256) * (5*sin(2*theta) + 12*sin(4*theta) + 33*sin(6*theta));
        elseif (m == 1); dPnm = sqrt(21)/128 * (5*cos(2*theta) + 24*cos(4*theta) + 99*cos(6*theta));
        elseif (m == 2); dPnm = sqrt(105/2)/256 * (-17*sin(2*theta) - 12*sin(4*theta) + 99*sin(6*theta));
        elseif (m == 3); dPnm = sqrt(105/2)*3/64 * sin(theta) .* (3*sin(3*theta) + 11*sin(5*theta));
        elseif (m == 4); dPnm = sqrt(7)*3/32 * (47*cos(theta) + 33*cos(3*theta)) .* sin(theta).^3;
        elseif (m == 5); dPnm = sqrt(77/2)*3/8 * (2 + 3*cos(2*theta)) .* sin(theta).^4;
        elseif (m == 6); dPnm = sqrt(231/2)*3/8 * cos(theta) .* sin(theta).^5;
        else;            dPnm = 0;
        end
    elseif (n == 7)
        if     (m == 0); dPnm = (-7/1024) * (25*sin(theta) + 81*sin(3*theta) + 165*sin(5*theta) + 429*sin(7*theta));
        elseif (m == 1); dPnm = sqrt(7)/2048 * (25*cos(theta) + 243*cos(3*theta) + 825*cos(5*theta) + 3003*cos(7*theta));
        elseif (m == 2); dPnm = -sqrt(21/2)/1024 * (75*sin(theta) + 171*sin(3*theta) + 55*sin(5*theta) - 1001*sin(7*theta));
        elseif (m == 3); dPnm = sqrt(21)/1024 * sin(theta) .* (45*sin(2*theta) + 396*sin(4*theta) + 1001*sin(6*theta));
        elseif (m == 4); dPnm = sqrt(231)/128 * (81 + 148*cos(2*theta) + 91*cos(4*theta)) .* sin(theta).^3;
        elseif (m == 5); dPnm = sqrt(231)/128 * (149*cos(theta) + 91*cos(3*theta)) .* sin(theta).^4;
        elseif (m == 6); dPnm = sqrt(3003/2)/32 * (5 + 7*cos(2*theta)) .* sin(theta).^5;
        elseif (m == 7); dPnm = sqrt(429)*7/32 * cos(theta) .* sin(theta).^6;
        else;            dPnm = 0;
        end
    elseif (n == 8)
        if     (m == 0); dPnm = (-9/2048) * (70*sin(2*theta) + 154*sin(4*theta) + 286*sin(6*theta) + 715*sin(8*theta));
        elseif (m == 1); dPnm = (3/1024) * (35*cos(2*theta) + 154*cos(4*theta) + 429*cos(6*theta) + 1430*cos(8*theta));
        elseif (m == 2); dPnm = sqrt(35/2)*3/512 * (-16*sin(2*theta) - 22*sin(4*theta) + 143*sin(8*theta));
        elseif (m == 3); dPnm = sqrt(1155)*3/512 * sin(theta) .* (3*sin(3*theta) + 13*(sin(5*theta) + 2*sin(7*theta)));
        elseif (m == 4); dPnm = sqrt(77)*3/128 * (138*cos(theta) + 117*cos(3*theta) + 65*cos(5*theta)) .* sin(theta).^3;
        elseif (m == 5); dPnm = sqrt(1001)*3/64 * (11 + 19*cos(2*theta) + 10*cos(4*theta)) .* sin(theta).^4;
        elseif (m == 6); dPnm = sqrt(429/2)*3/16 * (9*cos(theta) + 5*cos(3*theta)) .* sin(theta).^5;
        elseif (m == 7); dPnm = sqrt(715)*3/32 * (3 + 4*cos(2*theta)) .* sin(theta).^6;
        elseif (m == 8); dPnm = sqrt(715)*3/16 * cos(theta) .* sin(theta).^7;
        else;            dPnm = 0;
        end
    elseif (n == 9)
        if     (m == 0); dPnm = (-45/32768) * (98*sin(theta) + 11*(28*sin(3*theta) + 52*sin(5*theta) + 91*sin(7*theta) + 221*sin(9*theta)));
        elseif (m == 1); dPnm = sqrt(5)*3/32768 * (98*cos(theta) + 11*(84*cos(3*theta) + 260*cos(5*theta) + 637*cos(7*theta) + 1989*cos(9*theta)));
        elseif (m == 2); dPnm = sqrt(55/2)*3/8192 * (-98*sin(theta) - 252*sin(3*theta) - 260*sin(5*theta) + 91*sin(7*theta) + 1989*sin(9*theta));
        elseif (m == 3); dPnm = sqrt(1155/2)*3/8192 * sin(theta) .* (14*sin(2*theta) + 130*sin(4*theta) + 390*sin(6*theta) + 663*sin(8*theta));
        elseif (m == 4); dPnm = sqrt(5005)*3/2048 * (198 + 375*cos(2*theta) + 298*cos(4*theta) + 153*cos(6*theta)) .* sin(theta).^3;
        elseif (m == 5); dPnm = sqrt(143/2)*15/1024 * (418*cos(theta) + 325*cos(3*theta) + 153*cos(5*theta)) .* sin(theta).^4;
        elseif (m == 6); dPnm = sqrt(2145/2)*3/256 * (65 + 108*cos(2*theta) + 51*cos(4*theta)) .* sin(theta).^5;
        elseif (m == 7); dPnm = sqrt(715/2)*3/512 * (295*cos(theta) + 153*cos(3*theta)) .* sin(theta).^6;
        elseif (m == 8); dPnm = sqrt(12155)*3/256 * (7 + 9*cos(2*theta)) .* sin(theta).^7;
        elseif (m == 9); dPnm = sqrt(12155/2)*9/128 * cos(theta) .* sin(theta).^8;
        else;            dPnm = 0;
        end
    elseif (n == 10)
        if     (m == 0);  dPnm = (-55/65536) * (294*sin(2*theta) + 13*(48*sin(4*theta) + 81*sin(6*theta) + 136*sin(8*theta) + 323*sin(10*theta)));
        elseif (m == 1);  dPnm = sqrt(55)/32768 * (294*cos(2*theta) + 13*(96*cos(4*theta) + 243*cos(6*theta) + 544*cos(8*theta) + 1615*cos(10*theta)));
        elseif (m == 2);  dPnm = sqrt(165)/65536 * (-1666*sin(2*theta) + 13*(-208*sin(4*theta) - 171*sin(6*theta) + 136*sin(8*theta) + 1615*sin(10*theta)));
        elseif (m == 3);  dPnm = sqrt(2145/2)/8192 * sin(theta) .* (98*sin(3*theta) + 450*sin(5*theta) + 17*(63*sin(7*theta) + 95*sin(9*theta)));
        elseif (m == 4);  dPnm = sqrt(2145)/4096 * (4917*cos(theta) + 4455*cos(3*theta) + 3349*cos(5*theta) + 1615*cos(7*theta)) .* sin(theta).^3;
        elseif (m == 5);  dPnm = sqrt(429/2)*5/1024 * (572 + 1045*cos(2*theta) + 748*cos(4*theta) + 323*cos(6*theta)) .* sin(theta).^4;
        elseif (m == 6);  dPnm = sqrt(2145/2)/2048 * (5278*cos(theta) + 3859*cos(3*theta) + 1615*cos(5*theta)) .* sin(theta).^5;
        elseif (m == 7);  dPnm = sqrt(36465/2)/512 * (135 + 218*cos(2*theta) + 95*cos(4*theta)) .* sin(theta).^6;
        elseif (m == 8);  dPnm = sqrt(12155)/512 * (193*cos(theta) + 95*cos(3*theta)) .* sin(theta).^7;
        elseif (m == 9);  dPnm = sqrt(230945/2)/128 * (4 + 5*cos(2*theta)) .* sin(theta).^8;
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