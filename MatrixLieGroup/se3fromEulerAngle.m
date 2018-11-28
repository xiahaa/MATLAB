function se3 = se3fromEulerAngle(roll, pitch, yaw, x, y, z)
% construct se3 from euler angles and translation
%
rot = euler2rot(roll, pitch, yaw);
rotVec = rot2vec(rot);
se3 = [rotVec' x y z];
end
