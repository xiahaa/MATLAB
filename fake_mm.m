function fake_mm

clc;
clear all;
close all;

path = mfilename('fullpath');
i = findstr(path,'\');
folder = path(1:i(end));
cd(folder);

addpath('./MatrixLieGroup');
addpath('./3rdparty/rpnp1.0/code3/func')
addpath('./3rdparty/rpnp1.0/code3/epnp')
addpath('./3rdparty/rpnp1.0/code3/lhm')
addpath('./3rdparty/rpnp1.0/code3/sp')
addpath('./3rdparty/rpnp1.0/code3/sp')

addpath('./3rdparty/IPPE-master/IPPE-master/matlab')
addpath(genpath('./3rdparty/IPPE-master/IPPE-master/matlab'))

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
msize = 0.100;
se3_corners = generateMarkerCorners(se3_markers, msize);

%% fake camera position and traslation
numFrames = 30;
framePose = [];
i = 1;
dataCnt = 1;
isRandom = false;
if isRandom == false
    se3_fix = generateFixSE3();
    numFrames = size(se3_fix,1);
end
rerr = [];
err = [];
depth = [];
avgerr = [];
Tf = {};
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
   rot = vec2rot(cam_se3(1,1:3)');
   t = cam_se3(1,4:6);   
   
   rvs = [];
   tvs = [];
   depth = [depth;t(3)];
   
    t0 = se3_markers(1,4:6)';
    tc = Rc * rot * (t0 - t');
    r0 = euler2rot(se3_markers(1,1),se3_markers(1,2),se3_markers(1,3));
    rc = Rc*rot*r0';   
    Ts{i} = [rc tc;0 0 0 1];
   
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
        rotm = euler2rot(se3_markers(k,1),se3_markers(k,2),se3_markers(k,3));
        rvec = rot2vec(Rc*rot*rotm');
        
        %% test r, t
        rot_rec = vec2rot(rvec);
        pm = [-msize*0.5 msize*0.5 0; msize*0.5  msize*0.5 0; msize*0.5  -msize*0.5 0; -msize*0.5  -msize*0.5 0]';
        
        pdebug = rotm'*pm+repmat(tm,1,4);
        
        rvs = [rvs;rvec'];
        tvs = [tvs;tvec'];
        
        pc = rot_rec * pm + repmat(tvec,1,4);
        pc(:,1) = pc(:,1) ./ pc(3,1);pc(:,2) = pc(:,2) ./ pc(3,2);pc(:,3) = pc(:,3) ./ pc(3,3);pc(:,4) = pc(:,4) ./ pc(3,4);
        ppi = uint32(K * pc);
        plot(ppi(1,1), ppi(2,1),'w*');
        plot(ppi(1,2), ppi(2,2),'w*');
        plot(ppi(1,3), ppi(2,3),'w*');
        plot(ppi(1,4), ppi(2,4),'w*');        
   end   

   %% use rpnp as a pnp solver
   homogMethod = 'Harker';  %Harker DLT 
   opts.measureTiming = true;
   opts.withPoseRefinement = true; 
    frameerr = [];
    for j = 1:1:size(markersInView,1)
        pm = [-msize*0.5 msize*0.5 0; msize*0.5  msize*0.5 0; msize*0.5  -msize*0.5 0; -msize*0.5  -msize*0.5 0]';
        
        % 1st normarlize, 2*4
        imgpts = corners_total((j-1)*4+1:j*4,:)';
        
        imgpts(1,:) = (imgpts(1,:) - u0) ./ fx;
        imgpts(2,:) = (imgpts(2,:) - v0) ./ fy;
        
        %% TEST IPPE
%        [IPPEPoses, IPPERefinedPoses] = perspectiveIPPE(pm,imgpts,homogMethod,opts);
%         R1 = IPPEPoses.R1;
%         t1 = IPPEPoses.t1;
        
        [R1,t1]= RPnP(pm,imgpts);
%         imgpts = imgpts + randn(size(imgpts,1),size(imgpts,2));
%         [R1, t1] = lineRefinement(R1,t1,pm, imgpts, K)
%         tic
%         [R1,t1]= RPnP(pm,imgpts);
%         toc
%         rv_rpnp = rot2vec(R1);
        
        rerr = [rerr;norm(rot2vec(vec2rot(rvs(j,:)')'*R1))];
        frameerr = [frameerr;norm(t1 - tvs(j,:)')];
        err = [err;norm(t1 - tvs(j,:)')];
    end       
    avgerr = [avgerr;sum(sqrt(diag(frameerr*frameerr')))/size(frameerr,1)];

%    imshow(im,[]);
   pause(1);
   %% collect data
   fid = fopen(strcat('./data/',int2str(dataCnt),'.txt'),'wb');
   fid1 = fopen(strcat('./data/',int2str(dataCnt),'b.txt'),'w');
   fwrite(fid, size(markersInView,1), 'integer*4');
    for j = 1:1:size(markersInView,1)
        k = markersInView(j);
        fwrite(fid, k, 'integer*4');
        fprintf(fid1, '%d ', k);
        fwrite(fid, msize, 'single');
        fprintf(fid1, '%f ', msize);
        rvec = rvs(j,:);
        tvec = tvs(j,:);
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
    dataCnt = dataCnt + 1;
end

    mean(rerr)
    mean(err)

   fid = fopen(strcat('./data/','frameT','.txt'),'w');
    
   for i = 1:1:size(Ts,2)
       if ~isempty(Ts{i})
            fprintf(fid,'%f %f %f %f\n', Ts{i}(1,1),  Ts{i}(1,2),  Ts{i}(1,3),  Ts{i}(1,4));
            fprintf(fid,'%f %f %f %f\n', Ts{i}(2,1),  Ts{i}(2,2),  Ts{i}(2,3),  Ts{i}(2,4));
            fprintf(fid,'%f %f %f %f\n', Ts{i}(3,1),  Ts{i}(3,2),  Ts{i}(3,3),  Ts{i}(3,4));
            fprintf(fid,'%f %f %f %f\n', Ts{i}(4,1),  Ts{i}(4,2),  Ts{i}(4,3),  Ts{i}(4,4));
       end
   end
   fclose(fid);
end

function [Rrefine, trefine] = lineRefinement(R,t,objpts, imgpts, K)
    objptsEnh = [];
    lineId = [];
    lines = zeros(3,4);
    for i = 1:1:4
        if i+1>4
            i_1 = 1;
        else
            i_1 = i+1;
        end
        v1 = objpts(:,i_1) - objpts(:,i);
%         v1 = v1./norm(v1);
        s = linspace(0,1,10);
        for j = 1:1:length(s)
            objptsEnh = [objptsEnh objpts(:,i)+v1*s(j)];
            lineId = [lineId;i];
        end
        pt1 = [imgpts(:,i_1);1];
        pt2 = [imgpts(:,i);1];
        lines(:,i) = cross(pt1,pt2);
    end
    
    rV = rot2vec(R);
    x0 = [rV;t];
%     x0 = x0+randn(length(x0),1);
    
    Rinit = vec2rot(x0(1:3));
    tinit = x0(4:6);
    
    %% do the error computation and Jacobian computation
    OP = optimset;
    OP.Jacobian = 'on';
    %jacobians for STE are not currently implemented.
    OP.Display = 'on';
    %run levenberg marquardt:
    x = lsqnonlin(@(x)refineErr(x,objptsEnh,lines,lineId,K),x0,[],[],OP);
    
    Rrefine = vec2rot(x(1:3));
    trefine = x(4:6);
end

function [r,J] = refineErr(x,ps,lines,lineId,A)
rV = x(1:3);
t = x(4:6);
R = vec2rot(rV);

ps3D = R*ps;
ps3D(1,:) = ps3D(1,:)+t(1);
ps3D(2,:) = ps3D(2,:)+t(2);
ps3D(3,:) = ps3D(3,:)+t(3);

%% do projection
px = ps3D(1,:)./ps3D(3,:);
py = ps3D(2,:)./ps3D(3,:);

% px = A(1,1).*px + A(1,3);
% py = A(2,2).*py + A(2,3);
A1_1 = A(1,1);A1_2 = A(1,2);A1_3 = A(1,3);
A2_1 = A(2,1);A2_2 = A(2,2);A2_3 = A(2,3);
A3_1 = A(3,1);A3_2 = A(3,2);A3_3 = A(3,3);

%% compute line-to-distance err
r = zeros(length(lineId),1);
J = zeros(length(lineId),6);

rx = rV(1);
ry = rV(2);
rz = rV(3);
tx = t(1);
ty = t(2);
tz = t(3);

for i = 1:1:length(lineId)
    line = lines(:,lineId(i));
    err = [px(i) py(i) 1]*line;
    r(i) = err;
    
    l1 = line(1);
    l2 = line(2);
    l3 = line(3);

    Q1 = ps(1,i);
    Q2 = ps(2,i);
    Q3 = ps(3,i);
    
    J1 = l1*(Q2*((sin((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(rx)*conj(rz)*sign(rx))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(3/2) - (cos((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(rx)*conj(rz)*sign(rx))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2) - (conj(ry)*(cos((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2)) - 1))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2) + (2*abs(rx)*conj(rx)*conj(ry)*sign(rx)*(cos((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2)) - 1))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^2 + (sin((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(rx)*conj(rx)*conj(ry)*sign(rx))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(3/2)) - Q1*((2*conj(rx)*(cos((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2)) - 1))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2) + (sin((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(rx)*sign(rx))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2) - (sin((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(rx)*conj(rx)^2*sign(rx))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(3/2) - (2*abs(rx)*conj(rx)^2*sign(rx)*(cos((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2)) - 1))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^2) + Q3*((cos((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(rx)*conj(ry)*sign(rx))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2) - (conj(rz)*(cos((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2)) - 1))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2) - (sin((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(rx)*conj(ry)*sign(rx))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(3/2) + (2*abs(rx)*conj(rx)*conj(rz)*sign(rx)*(cos((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2)) - 1))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^2 + (sin((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(rx)*conj(rx)*conj(rz)*sign(rx))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(3/2))) + l2*(Q2*((sin((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(rx)*conj(ry)^2*sign(rx))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(3/2) - (sin((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(rx)*sign(rx))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2) + (2*abs(rx)*conj(ry)^2*sign(rx)*(cos((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2)) - 1))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^2) + Q3*((sin((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(rx)*conj(rx)*sign(rx))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(3/2) - (cos((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(rx)*conj(rx)*sign(rx))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2) - sin((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2) + (2*abs(rx)*conj(ry)*conj(rz)*sign(rx)*(cos((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2)) - 1))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^2 + (sin((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(rx)*conj(ry)*conj(rz)*sign(rx))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(3/2)) + Q1*((cos((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(rx)*conj(rz)*sign(rx))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2) - (conj(ry)*(cos((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2)) - 1))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2) - (sin((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(rx)*conj(rz)*sign(rx))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(3/2) + (2*abs(rx)*conj(rx)*conj(ry)*sign(rx)*(cos((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2)) - 1))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^2 + (sin((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(rx)*conj(rx)*conj(ry)*sign(rx))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(3/2))) + l3*(Q3*((sin((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(rx)*conj(rz)^2*sign(rx))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(3/2) - (sin((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(rx)*sign(rx))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2) + (2*abs(rx)*conj(rz)^2*sign(rx)*(cos((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2)) - 1))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^2) + Q2*(sin((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2) + (cos((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(rx)*conj(rx)*sign(rx))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2) - (sin((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(rx)*conj(rx)*sign(rx))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(3/2) + (2*abs(rx)*conj(ry)*conj(rz)*sign(rx)*(cos((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2)) - 1))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^2 + (sin((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(rx)*conj(ry)*conj(rz)*sign(rx))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(3/2)) + Q1*((sin((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(rx)*conj(ry)*sign(rx))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(3/2) - (cos((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(rx)*conj(ry)*sign(rx))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2) - (conj(rz)*(cos((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2)) - 1))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2) + (2*abs(rx)*conj(rx)*conj(rz)*sign(rx)*(cos((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2)) - 1))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^2 + (sin((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(rx)*conj(rx)*conj(rz)*sign(rx))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(3/2)));

    J2 = l2*(Q1*((cos((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(ry)*conj(rz)*sign(ry))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2) - (conj(rx)*(cos((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2)) - 1))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2) - (sin((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(ry)*conj(rz)*sign(ry))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(3/2) + (2*abs(ry)*conj(rx)*conj(ry)*sign(ry)*(cos((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2)) - 1))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^2 + (sin((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(ry)*conj(rx)*conj(ry)*sign(ry))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(3/2)) - Q2*((2*conj(ry)*(cos((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2)) - 1))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2) + (sin((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(ry)*sign(ry))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2) - (sin((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(ry)*conj(ry)^2*sign(ry))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(3/2) - (2*abs(ry)*conj(ry)^2*sign(ry)*(cos((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2)) - 1))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^2) + Q3*((sin((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(ry)*conj(rx)*sign(ry))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(3/2) - (cos((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(ry)*conj(rx)*sign(ry))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2) - (conj(rz)*(cos((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2)) - 1))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2) + (2*abs(ry)*conj(ry)*conj(rz)*sign(ry)*(cos((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2)) - 1))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^2 + (sin((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(ry)*conj(ry)*conj(rz)*sign(ry))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(3/2))) + l1*(Q1*((sin((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(ry)*conj(rx)^2*sign(ry))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(3/2) - (sin((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(ry)*sign(ry))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2) + (2*abs(ry)*conj(rx)^2*sign(ry)*(cos((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2)) - 1))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^2) + Q3*(sin((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2) + (cos((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(ry)*conj(ry)*sign(ry))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2) - (sin((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(ry)*conj(ry)*sign(ry))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(3/2) + (2*abs(ry)*conj(rx)*conj(rz)*sign(ry)*(cos((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2)) - 1))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^2 + (sin((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(ry)*conj(rx)*conj(rz)*sign(ry))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(3/2)) + Q2*((sin((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(ry)*conj(rz)*sign(ry))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(3/2) - (cos((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(ry)*conj(rz)*sign(ry))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2) - (conj(rx)*(cos((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2)) - 1))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2) + (2*abs(ry)*conj(rx)*conj(ry)*sign(ry)*(cos((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2)) - 1))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^2 + (sin((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(ry)*conj(rx)*conj(ry)*sign(ry))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(3/2))) + l3*(Q3*((sin((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(ry)*conj(rz)^2*sign(ry))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(3/2) - (sin((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(ry)*sign(ry))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2) + (2*abs(ry)*conj(rz)^2*sign(ry)*(cos((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2)) - 1))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^2) + Q1*((sin((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(ry)*conj(ry)*sign(ry))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(3/2) - (cos((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(ry)*conj(ry)*sign(ry))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2) - sin((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2) + (2*abs(ry)*conj(rx)*conj(rz)*sign(ry)*(cos((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2)) - 1))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^2 + (sin((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(ry)*conj(rx)*conj(rz)*sign(ry))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(3/2)) + Q2*((cos((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(ry)*conj(rx)*sign(ry))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2) - (conj(rz)*(cos((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2)) - 1))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2) - (sin((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(ry)*conj(rx)*sign(ry))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(3/2) + (2*abs(ry)*conj(ry)*conj(rz)*sign(ry)*(cos((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2)) - 1))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^2 + (sin((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(ry)*conj(ry)*conj(rz)*sign(ry))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(3/2)));

    J3 = l3*(Q1*((sin((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(rz)*conj(ry)*sign(rz))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(3/2) - (cos((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(rz)*conj(ry)*sign(rz))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2) - (conj(rx)*(cos((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2)) - 1))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2) + (2*abs(rz)*conj(rx)*conj(rz)*sign(rz)*(cos((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2)) - 1))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^2 + (sin((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(rz)*conj(rx)*conj(rz)*sign(rz))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(3/2)) - Q3*((2*conj(rz)*(cos((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2)) - 1))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2) + (sin((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(rz)*sign(rz))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2) - (sin((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(rz)*conj(rz)^2*sign(rz))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(3/2) - (2*abs(rz)*conj(rz)^2*sign(rz)*(cos((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2)) - 1))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^2) + Q2*((cos((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(rz)*conj(rx)*sign(rz))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2) - (conj(ry)*(cos((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2)) - 1))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2) - (sin((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(rz)*conj(rx)*sign(rz))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(3/2) + (2*abs(rz)*conj(ry)*conj(rz)*sign(rz)*(cos((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2)) - 1))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^2 + (sin((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(rz)*conj(ry)*conj(rz)*sign(rz))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(3/2))) + l1*(Q1*((sin((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(rz)*conj(rx)^2*sign(rz))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(3/2) - (sin((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(rz)*sign(rz))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2) + (2*abs(rz)*conj(rx)^2*sign(rz)*(cos((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2)) - 1))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^2) + Q2*((sin((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(rz)*conj(rz)*sign(rz))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(3/2) - (cos((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(rz)*conj(rz)*sign(rz))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2) - sin((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2) + (2*abs(rz)*conj(rx)*conj(ry)*sign(rz)*(cos((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2)) - 1))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^2 + (sin((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(rz)*conj(rx)*conj(ry)*sign(rz))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(3/2)) + Q3*((cos((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(rz)*conj(ry)*sign(rz))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2) - (conj(rx)*(cos((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2)) - 1))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2) - (sin((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(rz)*conj(ry)*sign(rz))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(3/2) + (2*abs(rz)*conj(rx)*conj(rz)*sign(rz)*(cos((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2)) - 1))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^2 + (sin((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(rz)*conj(rx)*conj(rz)*sign(rz))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(3/2))) + l2*(Q2*((sin((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(rz)*conj(ry)^2*sign(rz))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(3/2) - (sin((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(rz)*sign(rz))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2) + (2*abs(rz)*conj(ry)^2*sign(rz)*(cos((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2)) - 1))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^2) + Q1*(sin((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2) + (cos((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(rz)*conj(rz)*sign(rz))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2) - (sin((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(rz)*conj(rz)*sign(rz))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(3/2) + (2*abs(rz)*conj(rx)*conj(ry)*sign(rz)*(cos((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2)) - 1))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^2 + (sin((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(rz)*conj(rx)*conj(ry)*sign(rz))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(3/2)) + Q3*((sin((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(rz)*conj(rx)*sign(rz))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(3/2) - (cos((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(rz)*conj(rx)*sign(rz))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2) - (conj(ry)*(cos((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2)) - 1))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2) + (2*abs(rz)*conj(ry)*conj(rz)*sign(rz)*(cos((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2)) - 1))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^2 + (sin((abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(1/2))*abs(rz)*conj(ry)*conj(rz)*sign(rz))/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(3/2)));

    J4 = l1;

    J5 = l2;
    
    J6 = l3;
    
    J(i,:) = [J1 J2 J3 J4 J5 J6];
end

end



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
end

function se3_random = generateRandomSE3()
    roll = rand(1)*pi/4-pi/8;
    pitch = rand(1)*pi-pi/2;
    yaw = 0;
    se3_random(1:3) = euler2vec(roll,pitch,yaw)';
    se3_random(4) = (rand(1,1)*4)/10.0+0.6;
    se3_random(5) = rand(1,1)*4/10.0-0.2;
    se3_random(6) = rand(1,1)*3+2;
end

function crossangle = computeCrossAngle(se3_corners, se3_cam)
v1 = se3_corners(1,4:6) - se3_corners(2,4:6);
v1 = v1';
v2 = se3_corners(1,4:6) - se3_corners(4,4:6);
v2 = v2';
v3 = cross(v1,v2);
v3 = v3./norm(v3);

Rc = [1 0 0;0 -1 0;0 0 -1;];
rot_cam = vec2rot(se3_cam(1,1:3)');
v_z = Rc * rot_cam * [0 0 1]';
crossangle = acos(dot(v3,v_z));
end

function [corners_im] = fakeImage(corners, cam_se3, K)
rot = vec2rot(cam_se3(1,1:3)');
t = cam_se3(1,4:6);
corners_im = zeros(4,2);
Rc = [1 0 0;0 -1 0;0 0 -1;];
for i = 1:1:size(corners,1)
   p = K * Rc*(rot * (corners(i,4:6)' - t'));
   p = p./p(3);
%    p = K*p;
   corners_im(i,:) = [(p(1)) (p(2))];
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
end

%% verify is one marker is totally in current view
function inView = verifyMarkerInView(corners, cam_se3, K)
rot = vec2rot(cam_se3(1,1:3)');
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

end

function se3_corners = generateMarkerCorners(se3_markers, msize)
hsize = msize * 0.5;
corners = [-hsize hsize 0;
            hsize hsize 0;
            hsize -hsize 0;
            -hsize -hsize 0];

se3_corners = zeros(size(se3_markers,1)*4,6);        
        
%% do transfromation
for i = 1:1:size(se3_markers,1)
    rot = euler2rot(se3_markers(i,1),se3_markers(i,2),se3_markers(i,3));
    t = se3_markers(i,4:6);
    p = rot'*corners';
    p = p';
    p = p + repmat(t, 4, 1);
    vec = rot2vec(rot);
    se3_corners((i-1)*4+1:i*4,1:3) = repmat(vec',4,1);
    se3_corners((i-1)*4+1:i*4,4:6) = p;
end

end
