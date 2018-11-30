clc
clear all
close all

addpath('./solver/');
addpath('./data_gen/');
addpath('./evaluation/');
addpath('./stochastics/utils/');
addpath('../MatrixLieGroup/barfoot_tro14');
addpath('../quaternion');

%------------------------------------------------------
% opt == 1 : Use distribution of As and Bs
% opt == 2 : Use distribution of Xs
opt = 2;

%Editable Variables
%------------------------------------------------------
numBound = 100;	%Number of steps

breakPoint = [10 20 30 40 50 60 70 80 90 100];

gmean = [0;0;0;0;0;0];	%Gaussian Noise Mean

nstd = 0.01;  %Gaussian Noise standard deviation Range

shift = 0; %step shift

model = 1; %noise model

ElipseParam = [10, 10, 10];

trials = 1;

%%
index = 0; % index that records the number of pairs
dim = size(breakPoint,2);
N = 2; % Times of simulation
x = randn(6,1); X = expm(se3_vec(x));

for j = 1:N
    
    index = 0;
    
    for i = 1:numBound
        
        if    i == breakPoint(1) || i == breakPoint(2) || i == breakPoint(3)...
                || i == breakPoint(4) || i == breakPoint(5) || i == breakPoint(6)...
                || i == breakPoint(7) || i == breakPoint(8) || i == breakPoint(9)...
                || i == breakPoint(10)
            
            index = index + 1;
            numPair = i;
            
            %Trajectory Generation
            %------------------------------------------------------
            A = zeros(4,4,numPair);
            B = zeros(4,4,numPair);
            
            num = numPair + 2;
            t2 = (0:(2*pi)/((num + shift)):2*pi);
            twist = 0.0*sin(16*t2);
            
            trajParam = [.5, .5, .5, 0, 0];
            [A1, B1] = AB_genTraj(X, ElipseParam, trajParam, num/2, twist);
            
            trajParam = [.5, .5, .5, 0, 0.5*pi];
            [A2, B2] = AB_genTraj(X, ElipseParam, trajParam, num/2, twist);
            
            A = cat(3, A1, A2);
            B = cat(3, B1, B2);
            
            %-------------------------------------------------------
            if  opt == 1
                
                B = sensorNoise(B, gmean, nstd, 1);
                
            elseif opt ==2
                
                [A, X_dist] = A_NoiseX(B, X, gmean, nstd, 1);
                
            end
            
            for ii = 1:size(A,3)
                T1(ii,:,:) = B(:,:,ii);
                T2(ii,:,:) = A(:,:,ii);
            end
            
            %% Solution for X
            %-------------------------------------------------------
            numsam = size(T1,1);
            tic
            tsai(1:4,1:4,index) = sol_tsai_lenz(T1,T2,numsam);
            ttsai(index) = toc;
            
            tic
            tparkmartin(1:4,1:4,index) = sol_park_martin(T1,T2,numsam);
            ttparkmartin(index) = toc;
            
            tic
            thoraud(1:4,1:4,index) = sol_horaud(T1,T2,numsam);
            thoraud(index) = toc;
           
            tic
            thoraud_nlopt(1:4,1:4,index) = sol_horaud_nlopt(T1,T2,numsam);
            tthoraud_nlopt(index) = toc;
            
            tic
            tdq(1:4,1:4,index) = sol_dual_quaternion(T1,T2,numsam);
            tdq(index) = toc;
            
            tic
            tcvx(1:4,1:4,index) = sol_cvx1(T1,T2,numsam);
            ttcvx(index) = toc;
            
            tic
            tdphec(1:4,1:4,index) = sol_dphec(T1,T2,numsam);
            ttdphec(index) = toc;
            
            tic
            tidq(1:4,1:4,index) = sol_improved_dual_quaternion(T1,T2,numsam);
            ttidq(index) = toc;
            
            tic
            tata(1:4,1:4,index) = sol_adjoint_transformation_algo(T1,T2,numsam);
            ttata(index) = toc;
            
            tic
            tdual_cvx(1:4,1:4,index) = sol_dual_sdp_cvx(T1,T2,numsam);
            ttdual_cvx(index) = toc;
            
            tic
            tchou(1:4,1:4,index) = sol_chou(T1,T2,numsam);
            ttchou(index) = toc;
            
            % roterror and tranerror
            %--------------------------------------------------------
            % display('Kronecker')
            % error covariance matrix for the twist
            % the 6 parameters for the Lie Group
            %--------------------------------------------------------
            
            % trace of the error covariance matrix
            %-------------------------------------------------------
            rotError_tsai(index,j) = roterror(X, tsai(:,:,index));
            tranError_tsai(index,j) = tranerror(X, tsai(:,:,index));
            S_tsai = errorCovariance(A, B, tsai(:,:,index));
            tr_tsai(index,j) = trace(S_tsai);
            
            rotError_tparkmartin(index,j) = roterror(X, tparkmartin(:,:,index));
            tranError_tparkmartin(index,j) = tranerror(X, tparkmartin(:,:,index));
            S_tparkmartin = errorCovariance(A, B, tparkmartin(:,:,index));
            tr_tparkmartin(index,j) = trace(S_tparkmartin);
            
            rotError_thoraud(index,j) = roterror(X, thoraud(:,:,index));
            tranError_thoraud(index,j) = tranerror(X, thoraud(:,:,index));
            S_thoraud = errorCovariance(A, B, thoraud(:,:,index));
            tr_thoraud(index,j) = trace(S_thoraud);
            
            rotError_thoraud_nlopt(index,j) = roterror(X, thoraud_nlopt(:,:,index));
            tranError_thoraud_nlopt(index,j) = tranerror(X, thoraud_nlopt(:,:,index));
            S_thoraud_nlopt = errorCovariance(A, B, thoraud_nlopt(:,:,index));
            tr_thoraud_nlopt(index,j) = trace(S_thoraud_nlopt);
            
            rotError_tdq(index,j) = roterror(X, tdq(:,:,index));
            tranError_tdq(index,j) = tranerror(X, tdq(:,:,index));
            S_tdq = errorCovariance(A, B, tdq(:,:,index));
            tr_tdq(index,j) = trace(S_tdq);
            
            rotError_tcvx(index,j) = roterror(X, tcvx(:,:,index));
            tranError_tcvx(index,j) = tranerror(X, tcvx(:,:,index));
            S_tcvx = errorCovariance(A, B, tcvx(:,:,index));
            tr_tcvx(index,j) = trace(S_tcvx);
            
            rotError_tdphec(index,j) = roterror(X, tdphec(:,:,index));
            tranError_tdphec(index,j) = tranerror(X, tdphec(:,:,index));
            S_tdphec = errorCovariance(A, B, tdphec(:,:,index));
            tr_tdphec(index,j) = trace(S_tdphec);
            
            rotError_tidq(index,j) = roterror(X, tidq(:,:,index));
            tranError_tidq(index,j) = tranerror(X, tidq(:,:,index));
            S_tidq = errorCovariance(A, B, tidq(:,:,index));
            tr_tidq(index,j) = trace(S_tidq);
            
            rotError_tata(index,j) = roterror(X, tata(:,:,index));
            tranError_tata(index,j) = tranerror(X, tata(:,:,index));
            S_tata = errorCovariance(A, B, tata(:,:,index));
            tr_tata(index,j) = trace(S_tata);
            
            rotError_tdual_cvx(index,j) = roterror(X, tdual_cvx(:,:,index));
            tranError_tdual_cvx(index,j) = tranerror(X, tdual_cvx(:,:,index));
            S_tdual_cvx = errorCovariance(A, B, tdual_cvx(:,:,index));
            tr_tdual_cvx(index,j) = trace(S_tdual_cvx);
            
            rotError_tchou(index,j) = roterror(X, tchou(:,:,index));
            tranError_tchou(index,j) = tranerror(X, tchou(:,:,index));
            S_tchou = errorCovariance(A, B, tchou(:,:,index));
            tr_tchou(index,j) = trace(S_tchou);
        end
    end
end


%%
rotError = [sum(rotError_tsai,2)/N sum(rotError_tparkmartin,2)/N sum(rotError_thoraud,2)/N ...
            sum(rotError_thoraud_nlopt,2)/N sum(rotError_tdq,2)/N sum(rotError_tcvx,2)/N ...
            sum(rotError_tdphec,2)/N sum(rotError_tidq,2)/N sum(rotError_tata,2)/N ...
            sum(rotError_tdual_cvx,2)/N sum(rotError_tchou,2)/N];

tranError = [sum(tranError_tsai,2)/N ...
             sum(tranError_tparkmartin,2)/N ...
             sum(tranError_thoraud,2)/N ...
             sum(tranError_thoraud_nlopt,2)/N ...
             sum(tranError_tdq,2)/N ...
             sum(tranError_tcvx,2)/N ...
             sum(tranError_tdphec,2)/N ...
             sum(tranError_tidq,2)/N ...
             sum(tranError_tata,2)/N ...
             sum(tranError_tdual_cvx,2)/N ...
             sum(tranError_tchou,2)/N];

tr_Bar = [sum(tr_tsai,2)/N ...
            sum(tr_tparkmartin,2)/N ...
            sum(tr_thoraud,2)/N ...
            sum(tr_thoraud_nlopt,2)/N ...
            sum(tr_tdq,2)/N ...
            sum(tr_tcvx,2)/N ...
            sum(tr_tdphec,2)/N ...
            sum(tr_tidq,2)/N ...
            sum(tr_tdual_cvx,2)/N ...
            sum(tr_tchou,2)/N];


%%
cmap = jet(256);
figure
for i = 1:11
    plot([breakPoint],norm(rotError(:,i)),'Color',cmap(i*10,:));hold on;grid on;
end


xlabel('Number of AB Pairs')
ylabel('Rotation Error')
title('Rotation Error')

legend('Kron','Lie*','Quat','Dual Quat')


%%
figure
plot([breakPoint],tranErrorKronBar,'c')
hold on
plot([breakPoint],tranErrorLieBar,'r')
hold on
plot([breakPoint],tranErrorQuat,'b')
hold on
plot([breakPoint],tranErrorDual,'g')

xlabel('Number of AB Pairs')
ylabel('Rotation Error')
title('Translation Error')

legend('Kron','Lie*','Quat','Dual Quat')


%%
figure
plot([breakPoint],tr_kronBar,'c')
hold on
plot([breakPoint],tr_LieBar,'r')
hold on
plot([breakPoint],tr_quat,'b')
hold on
plot([breakPoint],tr_dual,'g')

xlabel('Number of AB Pairs')
ylabel('Trace')
title('Trace of Error Covariance')

legend('Kron','Lie*','Quat','Dual Quat')


%%