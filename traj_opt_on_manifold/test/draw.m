clc;close all;clear all;
[x,y,z]=sphere(30);
%     figure
h0=surf(x,y,z,'FaceAlpha',0.1,'FaceColor',[0.1,0.1,0.1],'EdgeColor','k');
%// Change your ranges here
minAzimuth = -180;
maxAzimuth = 180;
minElevation = 0;
maxElevation = 40;

%// Compute angles - assuming that you have already run the code for sphere
%// [x,y,z] = sphere;
%// x = x(11:end,:);
%// y = y(11:end,:); 
%// z = z(11:end,:);
theta = acosd(z);
phi = atan2d(y, x);

%%%%%// Begin highlighting logic
ind = (phi >= minAzimuth & phi <= maxAzimuth) & ...
      (theta >= minElevation & theta <= maxElevation); % // Find those indices
x2 = x; y2 = y; z2 = z; %// Make a copy of the sphere co-ordinates
x2(~ind) = NaN; y2(~ind) = NaN; z2(~ind) = NaN; %// Set those out of bounds to NaN

%%%%%// Draw our original sphere and then the region we want on top
r = 90;
surf(r.*x,r.*y,r.*z,'FaceColor','white','FaceAlpha',0.5); %// Base sphere
hold on;
surf(r.*x2,r.*y2,r.*z2,'FaceColor','red'); %// Highlighted portion
axis equal;
view(40,40); %// Adjust viewing angle for better view
axis off;
axis equal

%% sphere
figure
[x,y,z]=sphere(30);
r = 90;
surf(r.*x,r.*y,r.*z,'FaceColor','white','FaceAlpha',0.5); %// Base sphere