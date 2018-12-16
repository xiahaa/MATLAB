% main entry for comparing 4 SCF methods.
clc
clear all
close all

addpath('./solver/');
addpath('./data_gen/');
addpath('./evaluation/');
addpath('./stochastics/utils/');
addpath('../MatrixLieGroup/barfoot_tro14');
addpath('../quaternion');
addpath('./solver/atadq');


%Editable Variables
%------------------------------------------------------
% naaray = [10 20 30 40 50 60 70 80 90 100];
% naaray = [10 20 40 60 80 100];
naaray = [10 20 30 40];
gmean = [0;0;0;0;0;0];	%Gaussian Noise Mean
nstd = [0 0.01 0.05 0.1 0.2 0.5 0.8 1];  %Gaussian Noise standard deviation Range
nstd1 = [0 0.01 0.05 0.1 0.2 0.5 0.8 1];  %Gaussian Noise standard deviation Range
% nstd2 = [0 0.05 0.1 0.2 0.5 1];
nstd2 = [0.1 0.2 0.3 0.4 0.5 1];
sigma_r = [0 0.05 0.1 0.2 0.5 1];
sigma_t = [0 0.05 0.1 0.2 0.5 1];

usedstd = nstd2;

%%
N = 50; % Times of simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 advSolver1 = {
             @sol_cvx2_cmp};                             %% SCF

 usedsolver = advSolver1;
         
prefix1 = 'data/conv/convCmp';
prefix2 = 'data/SCF/scfCmp';
 
 
for j = 1:numel(usedstd)
    for i =  1:numel(naaray)
        numPair = naaray(i);
        %Trajectory Generation
        %------------------------------------------------------
        A = zeros(4,4,numPair);
        B = zeros(4,4,numPair); 
        num = numPair;
        
        noisylv = num2str(usedstd(j));
        noisylv = replace(noisylv,'.','_');
        filename1 = strcat(prefix1,'_',num2str(numPair), '_', noisylv, '.mat');
        dat = load(filename1);
        
        for k = 1:N
            X = dat.Xs{k};
            A = dat.As{k};
            B = dat.Bs{k};
            disp([i,j,k])
            %% re-arrange inputs
            for ii = 1:size(A,3)
                TB(ii,:,:) = B(:,:,ii);
                TA(ii,:,:) = A(:,:,ii);
            end
            %% Solution for X
            %-------------------------------------------------------
            numsam = size(TA,1);
            handle_sol = usedsolver{1};
            %% SCF
            [T1,T2,T3,T4,time1,time2,time3,time4] = handle_sol(TB,TA,numsam);
            xsol(1:4,1:4,1) = T1;
            xsol(1:4,1:4,2) = T2;
            xsol(1:4,1:4,3) = T3;
            xsol(1:4,1:4,4) = T4;
            flag(1) = 1;flag(2) = 1;flag(3) = 1;flag(4) = 1;
            
            xsols{k} = xsol;
            tsols{k} = [time1,time2,time3,time4];
            Xs{k} = X;
            flags{k} = flag;
        end
        noisylv = num2str(usedstd(j));
        noisylv = replace(noisylv,'.','_');
        filename = strcat(prefix2,'_',num2str(numPair), '_', noisylv, '.mat');
        save(filename,'xsols','tsols','Xs','flags');
    end
end
disp('finishing....');
