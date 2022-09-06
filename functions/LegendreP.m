% Returns the value of Schmidt-seminormalized Legendre function of degree n and order m given the angle theta
function Pnm = LegendreP(n, m, theta, UNNORM)  % n ( or l) = degree m = order
    if ~exist('UNNORM', 'var'); UNNORM = 0; end

    if (n == 1)
        if     (m == 0); Pnm = cos(theta);
        elseif (m == 1); Pnm = sin(theta);
        else;            Pnm = 0;
        end
    elseif (n == 2)
        if     (m == 0); Pnm = (1/4) * (1 + 3*cos(2*theta));
        elseif (m == 1); Pnm = (3) * cos(theta) .* sin(theta);
        elseif (m == 2); Pnm = (3) * sin(theta).^2;
        else;            Pnm = 0;
        end
    elseif (n == 3)
        if     (m == 0); Pnm = (1/8) * (3*cos(theta) + 5*cos(3*theta));
        elseif (m == 1); Pnm = (3/8) * (sin(theta) + 5*sin(3*theta));
        elseif (m == 2); Pnm = (-15) * cos(theta) .* sin(theta).^2;
        elseif (m == 3); Pnm = (15) * sin(theta).^3;
        else;            Pnm = 0;
        end
    elseif (n == 4)
        if     (m == 0); Pnm = (1/64) * (9 + 20*cos(2*theta) + 35*cos(4*theta));
        elseif (m == 1); Pnm = (5/16) * (2*sin(2*theta) + 7*sin(4*theta));
        elseif (m == 2); Pnm = (15/4) * (5 + 7*cos(2*theta)) .* sin(theta).^2;
        elseif (m == 3); Pnm = (105) * cos(theta) .* sin(theta).^3;
        elseif (m == 4); Pnm = (105) * sin(theta).^4;
        else;            Pnm = 0;
        end
    elseif (n == 5)
        if     (m == 0); Pnm = (1/128) * (30*cos(theta) + 35*cos(3*theta) + 63*cos(5*theta));
        elseif (m == 1); Pnm = (15/128) * (2*sin(theta) + 7*(sin(3*theta) + 3*sin(5*theta)));
        elseif (m == 2); Pnm = (105/8) * (5*cos(theta) + 3*cos(3*theta)) .* sin(theta).^2;
        elseif (m == 3); Pnm = (105/4) * (7 + 9*cos(2*theta)) .* sin(theta).^3;
        elseif (m == 4); Pnm = (945) * cos(theta) .* sin(theta).^4;
        elseif (m == 5); Pnm = (945) * sin(theta).^5;
        else;            Pnm = 0;
        end
    elseif (n == 6)
        if     (m == 0); Pnm = (1/512) * (50 + 105*cos(2*theta) + 126*cos(4*theta) + 231*cos(6*theta));
        elseif (m == 1); Pnm = (21/256) * (5*sin(2*theta) + 12*sin(4*theta) + 33*sin(6*theta));
        elseif (m == 2); Pnm = (105/64) * (35 + 60*cos(2*theta) + 33*cos(4*theta)) .* sin(theta).^2;
        elseif (m == 3); Pnm = (315/8) * (21*cos(theta) + 11*cos(3*theta)) .* sin(theta).^3;
        elseif (m == 4); Pnm = (945/4) * (9 + 11*cos(2*theta)) .* sin(theta).^4;
        elseif (m == 5); Pnm = (10395) * cos(theta) .* sin(theta).^5;
        elseif (m == 6); Pnm = (-10395) * sin(theta).^6;
        else;            Pnm = 0;
        end
    elseif (n == 7)
        if     (m == 0); Pnm = (1/1024) * (175*cos(theta) + 189*cos(3*theta) + 231*cos(5*theta) + 429*cos(7*theta));
        elseif (m == 1); Pnm = (7/1024) * (25*sin(theta) + 81*sin(3*theta) + 165*sin(5*theta) + 429*sin(7*theta));
        elseif (m == 2); Pnm = (63/128) * (350*cos(theta) + 275*cos(3*theta) + 143*cos(5*theta)) .* sin(theta).^2;
        elseif (m == 3); Pnm = (315/64) * (189 + 308*cos(2*theta) + 143*cos(4*theta)) .* sin(theta).^3;
        elseif (m == 4); Pnm = (3465/8) * (27*cos(theta) + 13*cos(3*theta)) .* sin(theta).^4;
        elseif (m == 5); Pnm = (10395/4) * (11 + 13*cos(2*theta)) .* sin(theta).^5;
        elseif (m == 6); Pnm = (-135135) * cos(theta) .* sin(theta).^6;
        elseif (m == 7); Pnm = (135135) * sin(theta).^7;
        else;            Pnm = 0;
        end
    elseif (n == 8)
        if     (m == 0); Pnm = (1/16384) * (1225 + 2520*cos(2*theta) + 2772*cos(4*theta) + 3432*cos(6*theta) + 6435*cos(8*theta));
        elseif (m == 1); Pnm = (9/2048) * (70*sin(2*theta) + 154*sin(4*theta) + 286*sin(6*theta) + 715*sin(8*theta));
        elseif (m == 2); Pnm = (315/2048) * (35 + 64*cos(2*theta) + 44*cos(4*theta) - 143*cos(8*theta));
        elseif (m == 3); Pnm = (3465/128) * (126*cos(theta) + 91*cos(3*theta) + 39*cos(5*theta)) .* sin(theta).^3;
        elseif (m == 4); Pnm = (10395/64) * (99 + 156*cos(2*theta) + 65*cos(4*theta)) .* sin(theta).^4;
        elseif (m == 5); Pnm = (135135/8) * (11*cos(theta) + 5*cos(3*theta)) .* sin(theta).^5;
        elseif (m == 6); Pnm = (-135135/4) * (13 + 15*cos(2*theta)) .* sin(theta).^6;
        elseif (m == 7); Pnm = (2027025) * cos(theta) .* sin(theta).^7;
        elseif (m == 8); Pnm = (2027025) * sin(theta).^8;
        else;            Pnm = 0;
        end
    elseif (n == 9)
        if     (m == 0); Pnm = (1/32768) * (4410*cos(theta) + 4620*cos(3*theta) + 143*(36*cos(5*theta) + 45*cos(7*theta) + 85*cos(9*theta)));
        elseif (m == 1); Pnm = (45/32768) * (98*sin(theta) + 11*(28*sin(3*theta) + 52*sin(5*theta) + 91*sin(7*theta) + 221*sin(9*theta)));
        elseif (m == 2); Pnm = (495/1024) * (735*cos(theta) + 637*cos(3*theta) + 455*cos(5*theta) + 221*cos(7*theta)) .* sin(theta).^2;
        elseif (m == 3); Pnm = (3465/512) * (462 + 819*cos(2*theta) + 546*cos(4*theta) + 221*cos(6*theta)) .* sin(theta).^3;
        elseif (m == 4); Pnm = (135135/128) * (66*cos(theta) + 45*cos(3*theta) + 17*cos(5*theta)) .* sin(theta).^4;
        elseif (m == 5); Pnm = (135135/64) * (143 + 220*cos(2*theta) + 85*cos(4*theta)) .* sin(theta).^5;
        elseif (m == 6); Pnm = (675675/8) * (39*cos(theta) + 17*cos(3*theta)) .* sin(theta).^6;
        elseif (m == 7); Pnm = (2027025/4) * (15 + 17*cos(2*theta)) .* sin(theta).^7;
        elseif (m == 8); Pnm = (34459425) * cos(theta) .* sin(theta).^8;
        elseif (m == 9); Pnm = (34459425) * sin(theta).^9;
        else;            Pnm = 0;
        end
    elseif (n == 10)
        if     (m == 0);  Pnm = (1/131072) * (7938 + 16170*cos(2*theta) + 17160*cos(4*theta) + 19305*cos(6*theta) + 24310*cos(8*theta) + 46189*cos(10*theta));
        elseif (m == 1);  Pnm = (55/65536) * (294*sin(2*theta) + 13*(48*sin(4*theta) + 81*sin(6*theta) + 136*sin(8*theta) + 323*sin(10*theta)));
        elseif (m == 2);  Pnm = (495/16384) * (8085 + 15288*cos(2*theta) + 12740*cos(4*theta) + 8840*cos(6*theta) + 4199*cos(8*theta)) .* sin(theta).^2;
        elseif (m == 3);  Pnm = (6435/1024) * (49*(33*cos(theta) + 27*cos(3*theta) + 17*cos(5*theta)) + 323*cos(7*theta)) .* sin(theta).^3;
        elseif (m == 4);  Pnm = (45045/512) * (858 + 1485*cos(2*theta) + 918*cos(4*theta) + 323*cos(6*theta)) .* sin(theta).^4;
        elseif (m == 5);  Pnm = (135135/128) * (1430*cos(theta) + 935*cos(3*theta) + 323*cos(5*theta)) .* sin(theta).^5;
        elseif (m == 6);  Pnm = (675675/64) * (585 + 884*cos(2*theta) + 323*cos(4*theta)) .* sin(theta).^6;
        elseif (m == 7);  Pnm = (11486475/8) * (45*cos(theta) + 19*cos(3*theta)) .* sin(theta).^7;
        elseif (m == 8);  Pnm = (34459425/4) * (17 + 19*cos(2*theta)) .* sin(theta).^8;
        elseif (m == 9);  Pnm = (654729075) * cos(theta) .* sin(theta).^9;
        elseif (m == 10); Pnm = (-654729075) * sin(theta).^10;
        else;             Pnm = 0;
        end
    else; Pnm = 0;
    end
    
    if ~UNNORM && m ~= 0
        SchmidtNormFac = sqrt(2*factorial(n-m)/factorial(n+m));
        Pnm = Pnm * SchmidtNormFac;
    end
end