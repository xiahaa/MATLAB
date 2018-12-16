% apply a certain hand eye calibration method and save the results.
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

ElipseParam = [10, 10, 10];

%%
index = 0; % index that records the number of pairs
N = 50; % Times of simulation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
convSolver = {@sol_tsai_lenz, ...                       %% TSAI
    @sol_park_martin, ...                               %% LIE
    @sol_horaud, ...                                    %% QSEP
    @sol_andreff, ...                                   %% KR
    @HandEye_DQ, ...                                    %% DQ
    @sol_improved_dual_quaternion, ...                  %% IDQ
    };

convSolver1 = {
    @HandEye_DQ, ...                                    %% DQ
    @sol_dual_quaternion,
    };

advSolver = {
             @sol_adjoint_transformation_algo, ...      %% ATA
             @sol_horaud_nlopt, ...                     %% NLOPT
             @sol_cvx1, ...                             %% SOCP
             @sol_dphec, ...                            %% GLOBAL_POLY
             @sol_dual_sdp_cvx, ...                     %% DUAL_SDP
             @sol_cvx2, ...                             %% SCF
             @sol_manifold_opt_SE3};                    %% SE3OPT
         
%  advSolver1 = {
%              @sol_horaud_nlopt, ...                     %% NLOPT
%              @sol_cvx1, ...                             %% SOCP
%              @sol_dphec, ...                            %% GLOBAL_POLY
%              @sol_dual_sdp_cvx, ...                     %% DUAL_SDP
%              @sol_cvx2, ...                             %% SCF
%              @sol_manifold_opt_SE3};                    %% SE3OPT

 advSolver1 = {
             @sol_horaud_nlopt, ...                     %% NLOPT
             };                    


 usedsolver = advSolver1;
         
prefix1 = 'data/conv/convCmp';
prefix2 = 'data/adv/nlopt/convCmp';
 
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
                T1(ii,:,:) = B(:,:,ii);
                T2(ii,:,:) = A(:,:,ii);
            end
            %% Solution for X
            %-------------------------------------------------------
            numsam = size(T1,1);
            for kk = 1:size(usedsolver, 2)
                handle_sol = usedsolver{kk};
                tic
                TX = handle_sol(T1,T2,numsam);
                tsol(kk) = toc;
                if isempty(TX) == false
                    xsol(1:4,1:4,kk) = TX;
                    flag(kk) = 1;
                else
                    flag(kk) = 0;
                end
            end
            xsols{k} = xsol;
            tsols{k} = tsol;
            Xs{k} = X;
            flags{k} = flag;
            As{k} = A;
            Bs{k} = B;
        end
        noisylv = num2str(usedstd(j));
        noisylv = replace(noisylv,'.','_');
        filename2 = strcat(prefix2,'_',num2str(numPair), '_', noisylv, '.mat');
        save(filename2,'xsols','tsols','Xs','flags','As','Bs');
    end
end
disp('finishing....');
