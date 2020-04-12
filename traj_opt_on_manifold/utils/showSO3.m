function showSO3(Rdata,Rreg)
    Rdata = reshape(Rdata,3,3,[]);
    Rreg = reshape(Rreg,3,3,[]);
    N1 = size(Rdata,3);
    N2 = size(Rreg,3);
    xid = zeros(3,N1);
    xi = zeros(3,N2);
    for i = 1:N1
        xid(:,i)=logSO3(Rdata(:,:,i));
    end
    for i = 1:N2
        xi(:,i) = logSO3(Rreg(:,:,i));
    end
    v = [1;1;1];%/sqrt(3);
    x1 = [xid(1,:).*v(1);xid(2,:).*v(2);xid(3,:)*v(3)];
    x2 = [xi(1,:).*v(1);xi(2,:).*v(2);xi(3,:)*v(3)];
    
    [x,y,z]=sphere(30);
%     figure
    h0=surf(x,y,z,'FaceAlpha',0.1,'FaceColor',[0.1,0.1,0.1],'EdgeColor','k')
%     shading interp
    hold on;
    cmap = lines(2);
    h1=plot3(x1(1,:),x1(2,:),x1(3,:),'o','MarkerSize',8,'Color',cmap(1,:),'MarkerFaceColor',cmap(1,:));
    h2=plot3(x2(1,:),x2(2,:),x2(3,:),'-','LineWidth',2,'Color',cmap(2,:));
%     title('Regression');
    legend([h1,h2],{'Input','Regression'});

    axis off;
    axis equal;
    
    set(gca,'xtick',[]);
    set(gca,'ytick',[]);
    set(gca,'ztick',[]);
    
    figure
    l = 0.2;
    step = 1;
    j = 1;
    for i = 1:10:N2
        draw_cube_plot([(j-1)*step,0,0],l,l,l,Rreg(:,:,i));hold on;
        odx = (j-1)*step;
        ody = 0;
        odz = 0;
        udx(1) = Rreg(1,1,i);
        vdx(1) = Rreg(2,1,i);
        wdx(1) = Rreg(3,1,i);
        udy(1) = Rreg(1,2,i);
        vdy(1) = Rreg(2,2,i);
        wdy(1) = Rreg(3,2,i);
        udz(1) = Rreg(1,3,i);
        vdz(1) = Rreg(2,3,i);
        wdz(1) = Rreg(3,3,i);
        quiver3(odx, ody, odz, udx, vdx, wdx, l,  'r', 'ShowArrowHead', 'off', 'MaxHeadSize', l, 'AutoScale', 'off','LineWidth',2);
        quiver3(odx, ody, odz, udy, vdy, wdy, l,  'g', 'ShowArrowHead', 'off', 'MaxHeadSize', l, 'AutoScale', 'off','LineWidth',2);
        quiver3(odx, ody, odz, udz, vdz, wdz, l,  'b', 'ShowArrowHead', 'off', 'MaxHeadSize', l, 'AutoScale', 'off','LineWidth',2);
        j = j + 1;
    end 
    view(3);
    axis equal;
end

function draw_cube_plot(origin,X,Y,Z,R)
% CUBE_PLOT plots a cube with dimension of X, Y, Z.
%
% INPUTS:
% origin = set origin point for the cube in the form of [x,y,z].
% X      = cube length along x direction.
% Y      = cube length along y direction.
% Z      = cube length along z direction.
% color  = STRING, the color patched for the cube.
%         List of colors
%         b blue
%         g green
%         r red
%         c cyan
%         m magenta
%         y yellow
%         k black
%         w white
% OUPUTS:
% Plot a figure in the form of cubics.
%
% EXAMPLES
% cube_plot(2,3,4,'red')
%

% ------------------------------Code Starts Here------------------------------ %
% Define the vertexes of the unit cubic

ver = [1 1 -1;
       -1 1 -1;
       -1 1 1;
       1 1 1;
       -1 -1 1;
       1 -1 1;
       1 -1 -1;
      -1 -1 -1];

%  Define the faces of the unit cubic
fac = [
    1 2 3 4;
    4 3 5 6;
    6 7 8 5;
    1 2 8 7;
    6 7 1 4;
    2 3 5 8
    ];
alphas = [0,0.6,0,0.6,0,0];
color = jet(6);
cube = [ver(:,1)*X,ver(:,2)*Y,ver(:,3)*Z];
% cube = [ver(:,1)*X+origin(1),ver(:,2)*Y+origin(2),ver(:,3)*Z+origin(3)];
cube = (R*cube' + origin')';
for i = 1:size(fac,1)
    patch('Faces',fac(i,:),'Vertices',cube,'FaceColor',color(i,:),'EdgeColor','k','FaceAlpha',alphas(i), ...
        'LineWidth',1.5);hold on;
end
hold off;
end
% ------------------------------Code Ends Here-------------------------------- %