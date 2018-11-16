function test_on_calibration_synchonization
    
    mpath = mfilename('fullpath');
    [folder, name, ext] = fileparts(mpath);
    cd(folder);
    
    %% fake some data;
    addpath('./MatrixLieGroup/barfoot_tro14');
    
    ti0 = 0;
    tend = 10;
    tis = 1/200.0;
    tcs = 1/10.0;
    ti = ti0:tis:tend;
    tc = ti0:tcs:tend;
    
    numi = numel(ti);
    numc = numel(tc);
    
    plotflag = 0;
    num = 2000;
    if 0
        for k = 1:num
            roll = (rand(numi,1).*180 - 90).*pi./180.0;
            pitch = (rand(numi,1).*180 - 90).*pi./180.0;
            yaw = (rand(numi,1).*360 - 180).*pi./180.0;

            dR = euler2rot((rand(1).*180 - 90).*pi./180.0, ...
                           (rand(1).*180 - 90).*pi./180.0, ...
                           (rand(1).*360 - 180).*pi./180.0);

            j = 1;
            Rimu = cell(numi,1);
            Rcam = cell(numc,1);
            for i = 1:numi        
                %% just for validation
                R1 = euler2rot(roll(i),pitch(i),yaw(i));
                Rimu{i} = R1';
                if j <= numc && ti(i) >= tc(j) 
                    R2 = dR*R1;
                    Rcam{j} = R2';
                    j = j + 1;
                end
            end
            trandshift = rand(1)*6 - 3;
            tcs = tc + trandshift;

            %% now we have the faked data
            tic
%             [roughalignedt] = timeSynchronizationPar(ti,Rimu,tcs,Rcam,-5,5,1);
%             [roughalignedt] = timeSynchronizationPar(ti,Rimu,tcs,Rcam,roughalignedt-3,roughalignedt+3,0.3);
%             [roughalignedt] = timeSynchronizationPar(ti,Rimu,tcs,Rcam,roughalignedt-1,roughalignedt+1,0.1);

%             [roughalignedt] = timeSynchronizationPar(ti,Rimu,tcs,Rcam,roughalignedt-0.3,roughalignedt+0.3,0.03);
%             [roughalignedt] = timeSynchronizationPar(ti,Rimu,tcs,Rcam,roughalignedt-1,roughalignedt+1,0.01);
%             toc
%             [roughalignedt] = timeSynchronizationPar(ti,Rimu,tcs,Rcam,trandshift-1,trandshift+1,0.01);
%             tic
            roughalignedt = trandshift;
            ttt = 0.5;
            [finealignedt,bestalignimu,bestidx,bestidy,aligntime] = timeSynchronizationPar(ti,Rimu,tcs,Rcam,trandshift-ttt,trandshift+ttt,0.005);
            toc
            if plotflag == 1
                close all;
                figure(1)
                plot(aligntime,bestalignimu,'g-o');hold on;grid on;
                plot(aligntime,bestaligncam,'r');hold on;
                title('Time synchronization');
                xlabel('t: (s)');
                ylabel('angle: (rad)');
                legend('IMU','Vision');
            end

            %% now assume the time is fully synchronized, start estimate the relative rotation
            % first compute the relative rotation
            numSample = numel(bestalignimu)-1;
            rRimu = cell(numSample,1);
            rRcam = cell(numSample,1);
            for i = 1:numSample
                R1 = Rimu{bestidx(i)};
                R2 = Rimu{bestidx(i+1)};
                dRimu = R2'*R1;

                R3 = Rcam{bestidy(i)};
                R4 = Rcam{bestidy(i+1)};
                dRcam = R4'*R3;

                rRimu{i} = dRimu;
                rRcam{i} = dRcam;
    %             dRsam{i} = inv(dRcam)*dR*dRimu*dR';
            end
            save(strcat('./data/calibration/data',num2str(k),'.mat'),'rRimu','rRcam','numSample', ...
                                                                     'dR','trandshift','roughalignedt','finealignedt');
            if (abs(trandshift-finealignedt)>=0.1)
                warning(num2str(abs(trandshift-finealignedt)));
            end
        end
    else
        t1 = zeros(num,1);
        t2 = zeros(num,1);
        t3 = zeros(num,1);
        error = zeros(num,1);
        for i = 1:num
            i
            load(strcat('./data/calibration/data',num2str(i),'.mat'));
            %% sample is ready, do rotation optimization
            % use lie algebra to derive the optimization problem
            dRopt = so3optimizationRANSAC(rRimu, rRcam, numSample);
            
            t1(i) = trandshift;
            t2(i) = roughalignedt;
            t3(i) = finealignedt;
            error(i) = norm(rot2vec(dRopt'*dR));
        end
        figure(1)
        plot(abs(t2-t1),'b-s');hold on;grid on;
        plot(abs(t3-t1),'r-o');
        title('Time');
        xlabel('Test id');ylabel('Error: (s)');
        legend('Rough search','Fine search');
        
        figure(2);
        plot(error,'LineWidth',3);hold on;grid on;
        xlabel('Test id');ylabel('Rotation Error: (rad)');

    end
    
    
end

function varargout = timeSynchronizationPar(tx,Rx,ty,Ry,tmin,tmax,tsetp)
    ts = tmin:tsetp:tmax;
    mts = numel(ts);
    m = size(Ry,1);
    bestalignx = [];
    bestaligny = [];
    berror1 = -1e8;
    berror2 = 1e8;
    bestshift = 0;
    aligntime = [];
    bestidx = [];
    bestidy = [];
    
    for k = 1:mts
        tshift = ts(k);
        txalign = tx + ts(k);
        ids = ones(m,1).*-1;
        for i = 1:m
            
            ids(i) = funmin(ty(i), txalign);
%             idf = -1;   
%             idf = funmin(ty(i), txalign);
%             if idf ~= -1
%                 ids(i) = idf;
%                 break;
%             end
        end
%         iii = (i + 1:1:m)';
%         if isempty(iii)
%             continue;
%         else
%             ids(iii) = 20.*(iii-i)+idf;
% %             ids(ids > (numel(tx) + 20)) = -1;
%             ids(ids > (numel(tx))) = numel(tx);
%         end
%         
        valididx = ids(ids ~= -1);
        valididy = 1:m;valididy(ids == -1) = [];valididy=valididy';
%         
%         if numel(valididx) < 10
%             continue;
%         end
        
        radx = zeros(numel(valididx)-1,1);
        rady = zeros(numel(valididx)-1,1);
        for i = 2:numel(valididx)
            R1 = Rx{valididx(i-1)};
            R2 = Rx{valididx(i)};
            R3 = Ry{valididy(i-1)};
            R4 = Ry{valididy(i)};
            dRx = R2'*R1;
            dRy = R4'*R3;
            radx(i-1) = norm(rot2vec(dRx));
            rady(i-1) = norm(rot2vec(dRy));
        end
        
        error1 = numel(radx);
        error2 = sum(abs(radx-rady)) / error1;
        if error1 > berror1 || (error1 == berror1 && berror2 >= error2)
            bestshift = tshift;
            berror1 = error1;berror2 = error2;
            bestalignx = [0;radx];
            bestidx = valididx;
            bestidy = valididy;
            aligntime = txalign(valididx);
        end
    end

    varargout{1} = bestshift;
    varargout{2} = bestalignx;
    varargout{4} = bestidy;
    varargout{3} = bestidx;
    varargout{5} = aligntime;   
end

function varargout = timeSynchronization(tx,Rx,ty,Ry,tmin,tmax,tsetp)
    bestalignx = [];
    bestaligny = [];
    berror1 = -1e8;
    berror2 = 1e8;
    bestshift = 0;
    aligntime = [];
    bestidx = [];
    bestidy = [];
    m = size(Ry,1);
    
    for tshift = tmin:tsetp:tmax
        txalign = tx + tshift;
        ids = zeros(m,1);
        for i = 1:m
            ids(i) = funmin(ty(i), txalign);
        end
        
        valididx = ids(ids ~= -1);
        valididy = 1:m;valididy(ids == -1) = [];valididy=valididy';
        
        if isempty(valididx)
            continue;
        end
        
        radx = zeros(numel(valididx)-1,1);
        rady = zeros(numel(valididy)-1,1);
        
        for i = 2:numel(valididx)
            R1 = Rx{valididx(i-1)};
            R2 = Rx{valididx(i)};
            R3 = Ry{valididy(i-1)};
            R4 = Ry{valididy(i)};
            dRx = R2'*R1;
            dRy = R4'*R3;
            radx(i-1) = norm(rot2vec(dRx));
            rady(i-1) = norm(rot2vec(dRy));
        end
        
        error1 = numel(radx);
        error2 = sum(abs(radx-rady)) / error1;
        if error1 > berror1 || (error1 == berror1 && berror2 >= error2)
            bestshift = tshift;
            berror1 = error1;berror2 = error2;
            bestalignx = [0;radx];
            bestaligny = [0;rady];
            bestidx = valididx;
            bestidy = valididy;
            aligntime = txalign(valididx);
        end
    end
    varargout{1} = bestshift;
    varargout{2} = bestalignx;
    varargout{3} = bestaligny;
    varargout{4} = bestidx;
    varargout{5} = bestidy;
    varargout{6} = aligntime;   
end

function id = funmin(y,x)
    v = x-y;
    [val,id]=min(abs(v));
    if abs(val) > 1.2
        id = -1;
    end
end