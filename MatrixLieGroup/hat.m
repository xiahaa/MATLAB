function vechat = hat(vec)

validateattributes(vec, {'double'},{'ncols',1})

if size(vec,1) == 3
    % so3
    vechat = [0 -vec(3) vec(2);vec(3) 0 -vec(1);-vec(2) vec(1) 0];
elseif size(vec,1) == 6
    % se3
    vecso3 = [0 -vec(6) vec(5); ...
              vec(6) 0 -vec(4); ... 
              -vec(5) vec(4) 0];
    vechat = [vecso3 vec(1:3);[0 0 0 0]];
else
    error('not yet implement hat opertation!\n');
end
end