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
naaray = [10 20 40 60 80 100];

nstd2 = [0 0.05 0.1 0.2 0.5 1];

usedstd = nstd2;

%%
N = 100; % Times of simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 advSolver1 = {
             @sol_manifold_opt_SE3_cmp};                             %% SCF

 usedsolver = advSolver1;
         
 prefix = 'data/SE3/Cmp';
         
for j = 3:numel(usedstd)
    for i =  6:numel(naaray)
        numPair = naaray(i);
        %Trajectory Generation
        %------------------------------------------------------
        A = zeros(4,4,numPair);
        B = zeros(4,4,numPair); 
        num = numPair;
        
        for k = 1:N
            disp([i,j,k])
            %%  rand X
            x = randn(6,1); X = expm(se3_vec(x));
            %% A,B
            sigma_r = usedstd(j);
            sigma_t = usedstd(j);
            for kx = 1:numPair       
                B(:,:,kx) = vec2tran(randn(6,1));
                A(:,:,kx) = X*B(:,:,kx)*X^-1;
                
                RBn = get_rotation_noise(sigma_r);
                tBn = randn(3,1)*sigma_t;
                B(1:3,1:3,kx) = RBn*B(1:3,1:3,kx);
                B(1:3,4,kx) = B(1:3,4,kx) + tBn;
                
                RAn = get_rotation_noise(sigma_r);
                tAn = randn(3,1)*sigma_t;
                A(1:3,1:3,kx) = RAn*A(1:3,1:3,kx);
                A(1:3,4,kx) = A(1:3,4,kx) + tAn;
            end            

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
            [T1,T2,T3,T4] = handle_sol(TB,TA,numsam);
            xsol(1:4,1:4,1) = T1;
            xsol(1:4,1:4,2) = T2;
            xsol(1:4,1:4,3) = T3;
            xsol(1:4,1:4,4) = T4;
%             xsol(1:4,1:4,5) = T5;
%             xsol(1:4,1:4,6) = T6;
            flag(1) = 1;flag(2) = 1;flag(3) = 1;flag(4) = 1;
%             flag(4) = 1;flag(5) = 1;flag(6) = 1;
            
            xsols{k} = xsol;
%             tsols{k} = [1,1,1,1,1,1];
            tsols{k} = [1,1,1,1];
            Xs{k} = X;
            flags{k} = flag;
        end
        noisylv = num2str(usedstd(j));
        noisylv = replace(noisylv,'.','_');
        filename = strcat(prefix,'_',num2str(numPair), '_', noisylv, '.mat');
        save(filename,'xsols','tsols','Xs','flags');
    end
end
disp('finishing....');
