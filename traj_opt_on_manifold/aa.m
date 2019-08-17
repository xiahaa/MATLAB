function aa
    x = bb(1,2,7,3,4)

end
function x = bb(x0,a,b,c,d)
    x = fzero(@cc,x0);
    function y = cc(x0)
        y = x0^3*a+b*x0^2+c*x0+d;
    end
end