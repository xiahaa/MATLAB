function [ C ] = vec2rot( phi )
% VEC2ROT Build a rotation matrix using the exponential map
%
% input: 
%   phi: 3x1 vector 
%
% output: 
%   C: 3x3 rotation matrix
%

validateattributes(phi,{'double'},{'size',[3,1]});

tolerance = 1e-12;
N = 10;

% Check for a small angle.
% 
angle = norm(phi);
if angle < tolerance
    % If the angle is small, fall back on the series representation.
    C = eye(3);
    xM = eye(3);
    cmPhi = hat(phi);
    for i = 1:N
        xM = xM * (cmPhi / i);
        C = C + xM;
    end
    %% normalize SO3
    C = normalizeSO3(C);
    %display('vec2rot.m:  used series method');
else
    axis = phi/angle;

    cp = cos(angle);
    sp = sin(angle);

    C = cp * eye(3) + (1 - cp) * axis * axis' + sp * hat(axis);
    %display('vec2rot.m:  used analytical method');
end

end

