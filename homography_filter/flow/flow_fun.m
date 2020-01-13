clc;clear; close all;
%% some parameters
fps = 30;
REF1 = 15;
REF2 = 7;
SAMP1 = 4;
SAMP2 = 3;
PSIZE = 16;
WIDTH = 320;
HEIGHT = 240;
XPATCHES = 10;
YPATCHES = 8;
XSPACING = 14;
YSPACING = 14;


XCENTRE = WIDTH/2;
YCENTRE = HEIGHT/2;
NUM_PATCHES = XPATCHES*YPATCHES; 

%% initsystem
goodMatch = ones(1,NUM_PATCHES);
imgpts0x = XCENTRE - XSPACING * (XPATCHES-1)/2 + ([1:XPATCHES]-1)*XSPACING;
imgpts0y = YCENTRE - YSPACING * (YPATCHES-1)/2 + ([1:YPATCHES]-1)*YSPACING;

%% todo
XMAX1 = XCENTRE - XSPACING*(XPATCHES-1)/2 - PSIZE*SAMP1/2 - REF1 - 3;    
YMAX1 = YCENTRE - YSPACING*(YPATCHES-1)/2 - PSIZE*SAMP1/2 - REF1 - 3;  
XMAX2 = XCENTRE - XSPACING*(XPATCHES-1)/2 - PSIZE*SAMP2/2 - REF2 - 3;    
YMAX2 = YCENTRE - YSPACING*(YPATCHES-1)/2 - PSIZE*SAMP2/2 - REF2 - 3;

Ric = eye(3);% assume no camera to imu rotation

%%
imgpts1x = imgpts0x; imgpts1y = imgpts0y; 
imgpts2x = imgpts0x; imgpts2y = imgpts0y; 
imgpts3x = imgpts0x; imgpts3y = imgpts0y; 
imgpts4x = imgpts0x; imgpts4y = imgpts0y; 

%% todo
for 1 
    %% image preFilter3 is actually the box filter 3x3
    img = loadimg();%% todo
    kernels = fspecial('average',3);
    preflt31 = imfilter(img, kernels);
    preflt15 = imfilter(img, kernels);
    kernelb = fspecial('average',32);
    kerne2b = fspecial('average',16);
    preflt31 = imfilter(preflt31, kernelb);
    preflt15 = imfilter(preflt15, kerne2b);
    
    %%
    imgpts1x = imgpts0x; imgpts1y = imgpts0y; 
    imgpts2x = imgpts0x; imgpts2y = imgpts0y; 
    imgpts3x = imgpts0x; imgpts3y = imgpts0y; 
    imgpts4x = imgpts0x; imgpts4y = imgpts0y; 
    
    I2A_OF(imgpre, imgcur, imgpts1x, imgpts1y, imgpts2x, imgpts2y, REF1, SAMP1);
    I2A_OF(imgpre, imgcur, imgpts1x, imgpts1y, imgpts2x, imgpts2y, SAMP2);
    
    imgpts4x = imgpts3x;
    imgpts4y = imgpts3y;
    copyPoints(imgpts3x, imgpts3y, imgpts4x, imgpts4y);
    I2A_OF(fltpre31, fltcur31, imgpts3x, imgpts3y, imgpts4x, imgpts4y, REF1, SAMP1);
    I2A_OF(fltpre15, fltcur15, imgpts3x, imgpts3y, imgpts4x, imgpts4y, REF2, SAMP2);
    
    I2A_OF(fltpre31, fltcur31, imgpts3x, imgpts3y, imgpts4x, imgpts4y, REF1, SAMP1);
    I2A_OF(fltpre15, fltcur15, imgpts3x, imgpts3y, imgpts4x, imgpts4y, REF2, SAMP2);
    
    %% one is warped of, one is no warped of
    calAveFlow(imgpts0x, imgpts0y, imgpts2x, imgpts2y, ave_ofwp);
    calAveFlow(imgpts3x, imgpts3y, imgpts4x, imgpts4y, ave_of);
end

function aveof = calAveFlow(p1x, p1y, p2x, p2y)
    aveof = zeros(1,3);
    aveof(1) = sum(sum(p2x - p1x));
    aveof(2) = sum(sum(p2y - p1y));
    signarray = ones(1,XPATCHES);
    signarray(1:round(signarray/2)) = -1;
    signarray = ones(YPATCHES,1)*YPATCHES;
    aveof(3) = sum(sum((p2y-p1y).*signarray));
    uxAve = sum(sum(p1x.*signarray));
    aveof(1) = aveof(1) / NUM_PATCHES;
    aveof(2) = aveof(2) / NUM_PATCHES;
    aveof(3) = aveof(3) / uxAve * FPS;
end

%%p1x, p1y: raw p2x, p2y: warpped by gyro step: 4 3 rshift: 15 7
function I2A_OF(imgpre, imgcur, p1x, p1y, p2x, p2y, rshift, step)
    PATCHLEN = step*PSIZE;
    rowmax1 = HEIGHT-PATCHLEN-rshift-1;
    colmax1 = WIDTH -PATCHLEN-rshift-1;
    scale = 100;
    for colp = 1:XPATCHES
        for rowp = 1:YPATCHES
            prex = (p2x(colp,rowp)-p1x(colp,rowp));
            prey = (p2y(colp,rowp)-p1y(colp,rowp));
            %% todo
            if (prex>100) prex = 100; end
            if (prex<-100) prex = -100; end
            if (prey>100) prey = 100; end
            if (prey<-100) prey = -100; end
            %% todo
            rowmax2 = HEIGHT-PATCHLEN-prey-1;
            colmax2 = WIDTH -PATCHLEN-prex-1;
            rowstart = int32(p1y(colp,rowp)-PATCHLEN/2);
            if (rowstart<rshift) rowstart = rshift; end
            if (rowstart<-prey) rowstart = -prey; end
            if (rowstart>rowmax1) rowstart = rowmax1; end
            if (rowstart>rowmax2) rowstart = rowmax2; end
            rowend = rowstart + PATCHLEN;
            colstart = int32(p1x[colp][rowp]-PATCHLEN/2);
            if (colstart<rshift) colstart = rshift; end
            if (colstart<-prex) colstart = -prex; end
            if (colstart>colmax1) colstart = colmax1; end
            if (colstart>colmax2) colstart = colmax2; end
            colend = colstart + PATCHLEN;
            A = 0;
            B = 0;
            C = 0;
            D = 0;
            E = 0;
            for row = rowstart:step:rowend
                for col = colstart:step:colend
                    f2_f1 = imgpre(col-rshift,row) - imgpre(col+rshift,row);
                    f4_f3 = imgpre(col,row-rshift) - imgpre(col,row+rshift);
                    f_f0  = imgcur(col+prex,row+prey) - imgpre(col,row);
                    A = A + f2_f1 * f2_f1 / scale;
                    B = B + f4_f3 * f2_f1 / scale;
                    C = C + f_f0  * f2_f1 * 2 / scale;
                    D = D + f4_f3 * f4_f3 / scale;
                    E = E + f_f0  * f4_f3 * 2 / scale;
                end
            end
            denom = (A*D) - (B*B);
            if denom ~= 0 
                dx = ((C*D)-(B*E)) / denom;
                dy = ((A*E)-(B*C)) / denom;
            end
            p2x(colp,rowp) = p1x(colp,rowp) + dx * rshift + prex;
            p2y(colp,rowp) = p1y(colp,rowp) + dy * rshift + prey;
        end
    end
end