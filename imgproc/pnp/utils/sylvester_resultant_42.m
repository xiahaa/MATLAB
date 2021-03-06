function out = sylvester_resultant_42(p1, p2)
% p1: x^4, x^3*y, x^2*y^2, x^2, x*y^3, x*y, y^4, y^2, 1
% p2: x^2 y^2 xy
    a1 = p1(1);
    a2 = p1(2);
    a3 = p1(3);
    a4 = p1(4);
    a5 = p1(5);
    a6 = p1(6);
    a7 = p1(7);
    a8 = p1(8);
    a9 = p1(9);

    b1 = p2(1);
    b2 = p2(2);
    b3 = p2(3);
    b4 = p2(4);

    %% [ x^8, x^6, x^4, x^2, 1]
    out = [ a1^2*b2^4 - a1*a2*b2^3*b3 - 2*a1*a3*b1*b2^3 + a1*a3*b2^2*b3^2 + 3*a1*a5*b1*b2^2*b3 - a1*a5*b2*b3^3 + 2*a1*a7*b1^2*b2^2 - 4*a1*a7*b1*b2*b3^2 + a1*a7*b3^4 + a2^2*b1*b2^3 - a2*a3*b1*b2^2*b3 - 2*a2*a5*b1^2*b2^2 + a2*a5*b1*b2*b3^2 + 3*a2*a7*b1^2*b2*b3 - a2*a7*b1*b3^3 + a3^2*b1^2*b2^2 - a3*a5*b1^2*b2*b3 - 2*a3*a7*b1^3*b2 + a3*a7*b1^2*b3^2 + a5^2*b1^3*b2 - a5*a7*b1^3*b3 + a7^2*b1^4, ...
            b4*a2^2*b2^3 - b4*a2*a3*b2^2*b3 - 4*b4*a2*a5*b1*b2^2 + b4*a2*a5*b2*b3^2 + 6*b4*a2*a7*b1*b2*b3 - b4*a2*a7*b3^3 + 2*a6*a2*b1*b2^3 - a8*a2*b1*b2^2*b3 - a4*a2*b2^3*b3 + 2*b4*a3^2*b1*b2^2 - 2*b4*a3*a5*b1*b2*b3 - 6*b4*a3*a7*b1^2*b2 + 2*b4*a3*a7*b1*b3^2 + 2*a8*a3*b1^2*b2^2 - 2*a4*a3*b1*b2^3 - a6*a3*b1*b2^2*b3 - 2*a1*b4*a3*b2^3 + a4*a3*b2^2*b3^2 + 3*b4*a5^2*b1^2*b2 - 3*b4*a5*a7*b1^2*b3 - ...
            2*a6*a5*b1^2*b2^2 - a8*a5*b1^2*b2*b3 + 3*a4*a5*b1*b2^2*b3 + a6*a5*b1*b2*b3^2 + 3*a1*b4*a5*b2^2*b3 - a4*a5*b2*b3^3 + 4*b4*a7^2*b1^3 - 2*a8*a7*b1^3*b2 + 2*a4*a7*b1^2*b2^2 + 3*a6*a7*b1^2*b2*b3 + a8*a7*b1^2*b3^2 + 4*a1*b4*a7*b1*b2^2 - 4*a4*a7*b1*b2*b3^2 - a6*a7*b1*b3^3 - 4*a1*b4*a7*b2*b3^2 + a4*a7*b3^4 - 2*a1*a8*b1*b2^3 + 2*a1*a4*b2^4 - a1*a6*b2^3*b3 + a1*a8*b2^2*b3^2, ...
            a3^2*b2^2*b4^2 - 2*a3*a4*b2^3*b4 - a3*a5*b2*b3*b4^2 - a3*a6*b2^2*b3*b4 - 6*a3*a7*b1*b2*b4^2 + a3*a7*b3^2*b4^2 + 4*a3*a8*b1*b2^2*b4 - 2*a9*a3*b1*b2^3 + a9*a3*b2^2*b3^2 + a4^2*b2^4 + 3*a4*a5*b2^2*b3*b4 - a4*a6*b2^3*b3 + 4*a4*a7*b1*b2^2*b4 - 4*a4*a7*b2*b3^2*b4 - 2*a4*a8*b1*b2^3 + a4*a8*b2^2*b3^2 + 3*a5^2*b1*b2*b4^2 - 4*a5*a6*b1*b2^2*b4 + a5*a6*b2*b3^2*b4 - 3*a5*a7*b1*b3*b4^2 - ...
            2*a5*a8*b1*b2*b3*b4 + 3*a9*a5*b1*b2^2*b3 - 2*a2*a5*b2^2*b4^2 - a9*a5*b2*b3^3 + a6^2*b1*b2^3 + 6*a6*a7*b1*b2*b3*b4 - a6*a7*b3^3*b4 - a6*a8*b1*b2^2*b3 + 2*a2*a6*b2^3*b4 + 6*a7^2*b1^2*b4^2 - 6*a7*a8*b1^2*b2*b4 + 2*a7*a8*b1*b3^2*b4 + 2*a9*a7*b1^2*b2^2 - 4*a9*a7*b1*b2*b3^2 + 2*a1*a7*b2^2*b4^2 + 3*a2*a7*b2*b3*b4^2 + a9*a7*b3^4 + a8^2*b1^2*b2^2 - 2*a1*a8*b2^3*b4 - a2*a8*b2^2*b3*b4 + ...
            2*a1*a9*b2^4 - a2*a9*b2^3*b3, a5^2*b2*b4^3 - 2*a5*a6*b2^2*b4^2 - a5*a7*b3*b4^3 - a5*a8*b2*b3*b4^2 + 3*a9*a5*b2^2*b3*b4 + a6^2*b2^3*b4 + 3*a6*a7*b2*b3*b4^2 - a6*a8*b2^2*b3*b4 - a9*a6*b2^3*b3 + 4*b1*a7^2*b4^3 - 6*b1*a7*a8*b2*b4^2 + a7*a8*b3^2*b4^2 + 2*a4*a7*b2^2*b4^2 + 4*a9*b1*a7*b2^2*b4 - 4*a9*a7*b2*b3^2*b4 - 2*a3*a7*b2*b4^3 + 2*b1*a8^2*b2^2*b4 - 2*a4*a8*b2^3*b4 - 2*a9*b1*a8*b2^3 ...
            + a9*a8*b2^2*b3^2 + 2*a3*a8*b2^2*b4^2 + 2*a4*a9*b2^4 - 2*a3*a9*b2^3*b4, ...
           (a9*b2^2 - a8*b2*b4 + a7*b4^2)^2];
 
end
