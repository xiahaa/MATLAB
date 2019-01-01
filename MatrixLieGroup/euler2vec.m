function vec = euler2vec(roll,pitch,yaw)
% from euler angle to rotation vector
% this is done by firstly construct a rotation matrix by ZYX order
% and then take the log 
rot = euler2rot(roll, pitch, yaw);
vec = rot2vec(rot);
end
