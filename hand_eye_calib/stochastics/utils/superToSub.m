function [ subPoses ] = superToSub( initialPose, supPoses )
%SUPERTOSUB Regenerates subscript poses from initial pose and relative
%superscript poses
%   Detailed explanation goes here
%   T^n_0 = T^n_(n-1)*T^(n-1)_(n-2)*...*T^1_(0) ----> sup pos
%   T^0_n = *T^0_(1)**T^1_(2)*...*T^(n-1)_(n)
%   equivalent to inverse a serise of transformation matrices
    n = size(supPoses,3);
    subPoses = zeros(4,4,n+1);
    subPoses(:,:,1) = initialPose;

    for i = 1:1:n
        subPoses(:,:,i+1) = subPoses(:,:,i)*supPoses(:,:,i);
    end

end

