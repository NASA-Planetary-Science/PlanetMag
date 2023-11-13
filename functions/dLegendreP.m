% Returns the value of the derivative of the Schmidt-seminormalized Legendre function of degree n 
% and order m given the angle theta
function dPnm = dLegendreP(n, m, theta, UNNORM)
    if ~exist('UNNORM', 'var'); UNNORM = 0; end
    
    if (n == 1)
        if     (m == 0); dPnm = -sin(theta);
        elseif (m == 1); dPnm = cos(theta);
        else;            dPnm = 0;
        end
    elseif (n == 2)
        if     (m == 0); dPnm = (-3) * cos(theta) .* sin(theta);
        elseif (m == 1); dPnm = (3) * cos(2*theta);
        elseif (m == 2); dPnm = (6) * cos(theta) .* sin(theta);
        else;            dPnm = 0;
        end
    elseif (n == 3)
        if     (m == 0); dPnm = -(3/8) * (sin(theta) + 5*sin(3*theta));
        elseif (m == 1); dPnm = (3/8) * (cos(theta) + 15*cos(3*theta));    
        elseif (m == 2); dPnm = (-15/4) * (sin(theta) - 3*sin(3*theta)); 
        elseif (m == 3); dPnm = (45) * cos(theta) .* sin(theta).^2;        
        else;            dPnm = 0;
        end
    elseif (n == 4)
        if     (m == 0); dPnm = (-5/16) * (2*sin(2*theta) + 7*sin(4*theta));
        elseif (m == 1); dPnm = (5/4) * (cos(2*theta) + 7*cos(4*theta));
        elseif (m == 2); dPnm = (15/4) * (-2*sin(2*theta) + 7*sin(4*theta));
        elseif (m == 3); dPnm = (105) * sin(theta) .* sin(3*theta);  
        elseif (m == 4); dPnm = (420) * cos(theta) .* sin(theta).^3;
        else;            dPnm = 0;
        end
    elseif (n == 5)
        if     (m == 0); dPnm = (-15/128) * (2*sin(theta) + 7*(sin(3*theta) + 3*sin(5*theta)));
        elseif (m == 1); dPnm = (15/128) * (2*cos(theta) + 21*(cos(3*theta) + 5*cos(5*theta)));
        elseif (m == 2); dPnm = (105/32) * (-2*sin(theta) - 3*sin(3*theta) + 15*sin(5*theta));
        elseif (m == 3); dPnm = (315/16) * sin(theta) .* (2*sin(2*theta) + 15*sin(4*theta));
        elseif (m == 4); dPnm = (945/2) * (3 + 5*cos(2*theta)) .* sin(theta).^3;
        elseif (m == 5); dPnm = (4725) * cos(theta) .* sin(theta).^4;
        else;            dPnm = 0;
        end
    elseif (n == 6)
        if     (m == 0); dPnm = (-21/256) * (5*sin(2*theta) + 12*sin(4*theta) + 33*sin(6*theta));
        elseif (m == 1); dPnm = (21/128) * (5*cos(2*theta) + 24*cos(4*theta) + 99*cos(6*theta));
        elseif (m == 2); dPnm = (105/128) * (-17*sin(2*theta) - 12*sin(4*theta) + 99*sin(6*theta));
        elseif (m == 3); dPnm = (945/16) * sin(theta) .* (3*sin(3*theta) + 11*sin(5*theta));
        elseif (m == 4); dPnm = (945/4) * (47*cos(theta) + 33*cos(3*theta)) .* sin(theta).^3;
        elseif (m == 5); dPnm = (10395) * (2 + 3*cos(2*theta)) .* sin(theta).^4;
        elseif (m == 6); dPnm = (62370) * cos(theta) .* sin(theta).^5;
        else;            dPnm = 0;
        end
    elseif (n == 7)
        if     (m == 0); dPnm = (-7/1024) * (25*sin(theta) + 81*sin(3*theta) + 165*sin(5*theta) + 429*sin(7*theta));
        elseif (m == 1); dPnm = (7/1024) * (25*cos(theta) + 243*cos(3*theta) + 825*cos(5*theta) + 3003*cos(7*theta));
        elseif (m == 2); dPnm = (-63/512) * (75*sin(theta) + 171*sin(3*theta) + 55*sin(5*theta) - 1001*sin(7*theta));
        elseif (m == 3); dPnm = (315/256) * sin(theta) .* (45*sin(2*theta) + 396*sin(4*theta) + 1001*sin(6*theta));
        elseif (m == 4); dPnm = (3465/16) * (81 + 148*cos(2*theta) + 91*cos(4*theta)) .* sin(theta).^3;
        elseif (m == 5); dPnm = (10395/8) * (149*cos(theta) + 91*cos(3*theta)) .* sin(theta).^4;
        elseif (m == 6); dPnm = (135135/2) * (5 + 7*cos(2*theta)) .* sin(theta).^5;
        elseif (m == 7); dPnm = (945945) * cos(theta) .* sin(theta).^6;
        else;            dPnm = 0;
        end
    elseif (n == 8)
        if     (m == 0); dPnm = (-9/2048) * (70*sin(2*theta) + 154*sin(4*theta) + 286*sin(6*theta) + 715*sin(8*theta));
        elseif (m == 1); dPnm = (9/512) * (35*cos(2*theta) + 154*cos(4*theta) + 429*cos(6*theta) + 1430*cos(8*theta));
        elseif (m == 2); dPnm = (315/256) * (-16*sin(2*theta) - 22*sin(4*theta) + 143*sin(8*theta));
        elseif (m == 3); dPnm = (10395/128) * sin(theta) .* (3*sin(3*theta) + 13*(sin(5*theta) + 2*sin(7*theta)));
        elseif (m == 4); dPnm = (10395/16) * (138*cos(theta) + 117*cos(3*theta) + 65*cos(5*theta)) .* sin(theta).^3;
        elseif (m == 5); dPnm = (135135/4) * (11 + 19*cos(2*theta) + 10*cos(4*theta)) .* sin(theta).^4;
        elseif (m == 6); dPnm = (405405) * (9*cos(theta) + 5*cos(3*theta)) .* sin(theta).^5;
        elseif (m == 7); dPnm = (2027025) * (3 + 4*cos(2*theta)) .* sin(theta).^6;
        elseif (m == 8); dPnm = (16216200) * cos(theta) .* sin(theta).^7;
        else;            dPnm = 0;
        end
    elseif (n == 9)
        if     (m == 0); dPnm = (-45/32768) * (98*sin(theta) + 11*(28*sin(3*theta) + 52*sin(5*theta) + 91*sin(7*theta) + 221*sin(9*theta)));
        elseif (m == 1); dPnm = (45/32768) * (98*cos(theta) + 11*(84*cos(3*theta) + 260*cos(5*theta) + 637*cos(7*theta) + 1989*cos(9*theta)));
        elseif (m == 2); dPnm = (495/4096) * (-98*sin(theta) - 252*sin(3*theta) - 260*sin(5*theta) + 91*sin(7*theta) + 1989*sin(9*theta));
        elseif (m == 3); dPnm = (10395/2048) * (14*sin(2*theta) + 130*sin(4*theta) + 390*sin(6*theta) + 663*sin(8*theta));
        elseif (m == 4); dPnm = (135135/256) * (198 + 375*cos(2*theta) + 298*cos(4*theta) + 153*cos(6*theta)) .* sin(theta).^3;
        elseif (m == 5); dPnm = (675675/128) * (418*cos(theta) + 325*cos(3*theta) + 153*cos(5*theta)) .* sin(theta).^4;
        elseif (m == 6); dPnm = (2027025/16) * (65 + 108*cos(2*theta) + 51*cos(4*theta)) .* sin(theta).^5;
        elseif (m == 7); dPnm = (2027025/8) * (295*cos(theta) + 153*cos(3*theta)) .* sin(theta).^6;
        elseif (m == 8); dPnm = (34459425/2) * (7 + 9*cos(2*theta)) .* sin(theta).^7;
        elseif (m == 9); dPnm = (310134825) * cos(theta) .* sin(theta).^8;
        else;            dPnm = 0;
        end
    elseif (n == 10)
        if     (m == 0);  dPnm = (-55/65536) * (294*sin(2*theta) + 13*(48*sin(4*theta) + 81*sin(6*theta) + 136*sin(8*theta) + 323*sin(10*theta)));
        elseif (m == 1);  dPnm = (55/32768) * (294*cos(2*theta) + 13*(96*cos(4*theta) + 243*cos(6*theta) + 544*cos(8*theta) + 1615*cos(10*theta)));
        elseif (m == 2);  dPnm = (495/32768) * (-1666*sin(2*theta) + 13*(-208*sin(4*theta) - 171*sin(6*theta) + 136*sin(8*theta) + 1615*sin(10*theta)));
        elseif (m == 3);  dPnm = (6435/2048) * sin(theta) .* (98*sin(3*theta) + 450*sin(5*theta) + 17*(63*sin(7*theta) + 95*sin(9*theta)));
        elseif (m == 4);  dPnm = (45045/512) * (4917*cos(theta) + 4455*cos(3*theta) + 3349*cos(5*theta) + 1615*cos(7*theta)) .* sin(theta).^3;
        elseif (m == 5);  dPnm = (675675/128) * (572 + 1045*cos(2*theta) + 748*cos(4*theta) + 323*cos(6*theta)) .* sin(theta).^4;
        elseif (m == 6);  dPnm = (675675/64) * (5278*cos(theta) + 3859*cos(3*theta) + 1615*cos(5*theta)) .* sin(theta).^5;
        elseif (m == 7);  dPnm = (11486475/8) * (135 + 218*cos(2*theta) + 95*cos(4*theta)) .* sin(theta).^6;
        elseif (m == 8);  dPnm = (34459425/4) * (193*cos(theta) + 95*cos(3*theta)) .* sin(theta).^7;
        elseif (m == 9);  dPnm = (654729075) * (4 + 5*cos(2*theta)) .* sin(theta).^8;
        elseif (m == 10); dPnm = (6547290750) * cos(theta) .* sin(theta).^9;
        else;             dPnm = 0;
        end
    else; dPnm = 0;
    end
    
    if ~UNNORM && m ~= 0
        SchmidtNormFac = sqrt(2*factorial(n-m)/factorial(n+m));
        dPnm = dPnm * SchmidtNormFac;
    end
end