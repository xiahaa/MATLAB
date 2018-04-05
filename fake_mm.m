function fake_mm
fx = 4000;
fy = 4000;
u0 = 1200;
v0 = 1000;
%% faked camera
K = [fx 0 u0;0 fy v0;0 0 1];
%% faked markers in 3D space
numMarkers = 7;
se3_markers = zeros(numMarkers,6);
se3_markers(1,:) = [0,  0,  0,  0,      0,      0];%% origo
se3_markers(2,:) = [0,  0,  0,  0.4,    -0.4,   0];%% 1
se3_markers(3,:) = [0,  0,  deg2rad(5),0.4,0.4,0];%% 2
se3_markers(4,:) = [0,  deg2rad(45),0,-0.5,0.1,1];%% 3
se3_markers(5,:) = [0,deg2rad(45),0,-0.5,-0.2,1.2];%% 4
se3_markers(6,:) = [0,deg2rad(-45),0,0.8,0.15,1];%% 5
se3_markers(7,:) = [0,deg2rad(-45),0,0.8,-0.15,1.1];%% 6

%% generate markers' corners in 3D space
msize = 0.155;
se3_corners = generateMarkerCorners(se3_markers, msize);

%% fake camera position and traslation
numFrames = 30;
framePose = [];
i = 1;
isRandom = true;
if isRandom == false
    se3_fix = generateFixSE3();
    numFrames = size(se3_fix,1);
end
while( i < numFrames)
    if isRandom == true
        se3_random = generateRandomSE3();
    else
        se3_random = se3_fix(i,:);
    end
    cam_se3 = se3_random;
    %% check how many markers are in view
    markersInView = [];
    for j = 1:1:numMarkers
        inView = verifyMarkerInView(se3_corners((j-1)*4+1:j*4,:), cam_se3, K);
        if inView == 1
            crossangle = computeCrossAngle(se3_corners((j-1)*4+1:j*4,:), cam_se3);
            if (crossangle < deg2rad(45))
                markersInView = [markersInView;j];
            end
        end
    end
    if size(markersInView,1) < 2 
        if isRandom == false
            i = i+1;
        end
        continue;
    end
    %% display an image
   im = uint8(zeros(v0*2,u0*2));
   imshow(im,[]);
   hold on;
   corners_total = [];
   
   Rc = [1 0 0;0 -1 0;0 0 -1;];
   rot = rodrigues(cam_se3(1,1:3));
   t = cam_se3(1,4:6);   
   
   for j = 1:1:size(markersInView,1)
        k = markersInView(j);
        [corners] = fakeImage(se3_corners((k-1)*4+1:k*4,:), cam_se3, K);
        
        plot([corners(1,1) corners(2,1)], [corners(1,2) corners(2,2)], '-');
        plot([corners(2,1) corners(3,1)], [corners(2,2) corners(3,2)], '-');
        plot([corners(3,1) corners(4,1)], [corners(3,2) corners(4,2)], '-');
        plot([corners(4,1) corners(1,1)], [corners(4,2) corners(1,2)], '-');        
%         imshow(im,[]);
        corners_total = [corners_total;corners];
        
        tm = se3_markers(k,4:6)';
        tvec = Rc * rot * (tm - t');
        rotm = eular2rot(se3_markers(k,1),se3_markers(k,2),se3_markers(k,3));
        rvec = rodrigues(Rc*rot*rotm');
        
        %% test r, t
        rot_rec = rodrigues(rvec);
        pm = [-msize*0.5 msize*0.5 0; msize*0.5  msize*0.5 0; msize*0.5  -msize*0.5 0; -msize*0.5  -msize*0.5 0]';
        
        pdebug = rotm'*pm+repmat(tm,1,4);
        
        pc = rot_rec * pm + repmat(tvec,1,4);
        pc(:,1) = pc(:,1) ./ pc(3,1);pc(:,2) = pc(:,2) ./ pc(3,2);pc(:,3) = pc(:,3) ./ pc(3,3);pc(:,4) = pc(:,4) ./ pc(3,4);
        ppi = uint32(K * pc);
        plot(ppi(1,1), ppi(2,1),'w*');
        plot(ppi(1,2), ppi(2,2),'w*');
        plot(ppi(1,3), ppi(2,3),'w*');
        plot(ppi(1,4), ppi(2,4),'w*');
        
   end   
%    imshow(im,[]);
   pause(1);
   %% collect data
   fid = fopen(strcat(int2str(i),'.txt'),'wb');
   fid1 = fopen(strcat(int2str(i),'b.txt'),'w');
   fwrite(fid, size(markersInView,1), 'integer*4');
    for j = 1:1:size(markersInView,1)
        k = markersInView(j);
        fwrite(fid, k, 'integer*4');
        fprintf(fid1, '%d ', k);
        fwrite(fid, msize, 'single');
        fprintf(fid1, '%f ', msize);

        fwrite(fid, rvec, 'single');
        fprintf(fid1, '%f ', rvec);
        fwrite(fid, tvec, 'single');
        fprintf(fid1, '%f ', tvec);
        fwrite(fid, 4, 'integer*4');
        fprintf(fid1, '%d ', 4);
        imgpts = corners_total((j-1)*4+1:j*4,:);
        fwrite(fid, imgpts(1,1), 'single'); fwrite(fid, imgpts(1,2), 'single');
        fprintf(fid1, '%f ', imgpts(1,:));
        fwrite(fid, imgpts(2,1), 'single'); fwrite(fid, imgpts(2,2), 'single');        fprintf(fid1, '%f ', imgpts(2,:));
        fwrite(fid, imgpts(3,1), 'single'); fwrite(fid, imgpts(3,2), 'single');        fprintf(fid1, '%f ', imgpts(3,:));
        fwrite(fid, imgpts(4,1), 'single'); fwrite(fid, imgpts(4,2), 'single');        fprintf(fid1, '%f ', imgpts(4,:));fprintf(fid1,'\n');

    end
    fclose(fid);
    fclose(fid1);
    i = i+1;
end
return

function se3_fix = generateFixSE3()
i = 1;
se3_fix = [];
se3 = zeros(1,6);
for range = 5:-1:1
    se3(1,4:6) = [0.5,0,range];
    yaw = 0;
    for pitch = -pi/2:pi/8:pi/2
        for roll = -pi/8:pi/16:pi/8
            se3(1,1:3) = euler2vec(roll,pitch,yaw)';
            se3_fix = [se3_fix;se3];
        end
    end
end
return

function se3_random = generateRandomSE3()
    roll = rand(1)*pi/4-pi/8;
    pitch = rand(1)*pi-pi/2;
    yaw = 0;
    se3_random(1:3) = euler2vec(roll,pitch,yaw)';
    se3_random(4) = (rand(1,1)*4)/10.0+0.6;
    se3_random(5) = rand(1,1)*4/10.0-0.2;
    se3_random(6) = rand(1,1)*3+2;
return

function crossangle = computeCrossAngle(se3_corners, se3_cam)
v1 = se3_corners(1,4:6) - se3_corners(2,4:6);
v1 = v1';
v2 = se3_corners(1,4:6) - se3_corners(4,4:6);
v2 = v2';
v3 = cross(v1,v2);
v3 = v3./norm(v3);

Rc = [1 0 0;0 -1 0;0 0 -1;];
rot_cam = rodrigues(se3_cam(1,1:3));
v_z = Rc * rot_cam * [0 0 1]';
crossangle = acos(dot(v3,v_z));
return

function [corners_im] = fakeImage(corners, cam_se3, K)
rot = rodrigues(cam_se3(1,1:3));
t = cam_se3(1,4:6);
corners_im = zeros(4,2);
Rc = [1 0 0;0 -1 0;0 0 -1;];
for i = 1:1:size(corners,1)
   p = K * Rc*(rot * (corners(i,4:6)' - t'));
   p = p./p(3);
%    p = K*p;
   corners_im(i,:) = [uint32(p(1)) uint32(p(2))];
end
% im(corners_im(1,2):corners_im(2,2),corners_im(1,1):corners_im(2,1)) = 255;
% im(corners_im(2,2):corners_im(3,2),corners_im(2,1):corners_im(3,1)) = 255;
% im(corners_im(3,2):corners_im(4,2),corners_im(3,1):corners_im(4,1)) = 255;
% im(corners_im(4,2):corners_im(1,2),corners_im(4,1):corners_im(1,1)) = 255;
% minu = min(corners_im(:,1));
% minv = min(corners_im(:,2));
% maxu = max(corners_im(:,1));
% maxv = max(corners_im(:,2));
% width = maxu - minu;
% height = maxv - minv;
% im(minv:1:maxv,minu:1:maxu) = 255;
return

%% verify is one marker is totally in current view
function inView = verifyMarkerInView(corners, cam_se3, K)
rot = rodrigues(cam_se3(1,1:3));
t = cam_se3(1,4:6);
inView = 1;
Rc = [1 0 0;0 -1 0;0 0 -1;];
margin = 20;
for i = 1:1:size(corners,1)
   p = Rc*(rot * (corners(i,4:6)' - t'));
   if (p(3) < 0)
       inView = 0;
       return
   end
   p = p./p(3);
   p = K*p;
  
   if p(1)<margin || p(1) > (K(1,3)*2-margin) || p(2) < margin || p(2) > (K(2,3)*2-margin)
       inView = 0;
       return;
   end
end

return

function se3_corners = generateMarkerCorners(se3_markers, msize)
hsize = msize * 0.5;
corners = [-hsize hsize 0;
            hsize hsize 0;
            hsize -hsize 0;
            -hsize -hsize 0];

se3_corners = zeros(size(se3_markers,1)*4,6);        
        
%% do transfromation
for i = 1:1:size(se3_markers,1)
    rot = eular2rot(se3_markers(i,1),se3_markers(i,2),se3_markers(i,3));
    t = se3_markers(i,4:6);
    p = rot'*corners';
    p = p';
    p = p + repmat(t, 4, 1);
    vec = rodrigues(rot);
    se3_corners((i-1)*4+1:i*4,1:3) = repmat(vec',4,1);
    se3_corners((i-1)*4+1:i*4,4:6) = p;
end

return

%% from degree to rad, OK
function rad = deg2rad(deg)
rad = deg * pi / 180.0;
return

%% from eular angle, translation to se3, OK
function se3 = fromEularAngle(roll, pitch, yaw, x, y, z)
rot = eular2rot(roll, pitch, yaw);
rotVec = rodrigues(rot);
se3 = [rotVec' x y z];
return

%% OK
function vec = euler2vec(roll,pitch,yaw)
rot = eular2rot(roll, pitch, yaw);
vec = rodrigues(rot);
return

%% from eular angle to rotation matrix, OK
function rot = eular2rot(roll, pitch, yaw)
rot1 = [cos(yaw) sin(yaw) 0;-sin(yaw) cos(yaw) 0; 0 0 1];
rot2 = [cos(pitch) 0 -sin(pitch);0 1 0; sin(pitch) 0 cos(pitch)];
rot3 = [1 0 0;0 cos(roll) sin(roll);0 -sin(roll) cos(roll)];

rot = rot3*rot2*rot1;

return

%% from rotation matrix to rotation vector, OK
function rotVec = rot2vec(rot)
theta = acos((trace(rot) - 1)*0.5);
[U S V] = svd(rot - eye(3,3));
axis = V(:,end);
rotVec = theta.*axis;
return

%% OK
function rot = vec2rot(vec)
theta = norm(vec);
vec = vec ./ (theta+1e-7);
rot = eye(3,3) + sin(theta) * skewsymm(vec) + (1-cos(theta))*skewsymm(vec)*skewsymm(vec);

% if theta < 1e-6
%     rot = eye(3,3);
%     return;
% end
% axis = vec ./ theta;
% rot = cos(theta).*eye(3,3)+(1-cos(theta)).*(axis*axis')+sin(theta)*skewsymm(axis);

return

