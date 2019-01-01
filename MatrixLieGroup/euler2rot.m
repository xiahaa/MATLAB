function rot = euler2rot(roll, pitch, yaw)
% from eular angle to rotation matrix, OK

rot1 = [cos(yaw) sin(yaw) 0;-sin(yaw) cos(yaw) 0; 0 0 1];
rot2 = [cos(pitch) 0 -sin(pitch);0 1 0; sin(pitch) 0 cos(pitch)];
rot3 = [1 0 0;0 cos(roll) sin(roll);0 -sin(roll) cos(roll)];

rot = rot3*rot2*rot1;

end
