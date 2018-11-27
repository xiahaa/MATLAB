function [ X, MeanA, MeanB, SigA, SigB, t_error ] = batchSolveNew(A, B, opt)
    %% Mixed version for solving AX = XB
    %%
    [a1,a2,a3]  = size(A);

    n_search = int16(2*10^2);

    if opt == 1
        [MeanA, ~] = mean_Taylor_1st( A ); %_mex
        [MeanB, ~] = mean_Taylor_1st( B ); %_mex
    elseif opt == 2
        MeanA = mean_Taylor_2nd( A, 0, n_search ); %_mex
        MeanB = mean_Taylor_2nd( B, 0, n_search ); %_mex
    elseif opt == 3
        [ MeanA, ~ ] = distibutionPropsMex( A );
        [ MeanB, ~ ] = distibutionPropsMex( B );
    elseif opt == 4
        MeanA = mean_Taylor_2nd_adv_recursive( A, 1, n_search ); 
        MeanB = mean_Taylor_2nd_adv_recursive( B, 1, n_search ); 
    end
    
    SigA = zeros(6,6);
    SigB = zeros(6,6);
    SigA = cov_SE3(MeanA, A, 1);
    SigB = cov_SE3(MeanB, B, 1);

%     for i = 1:size(A,3)
%         SigA = SigA + se3_vec(logm(MeanA^(-1)*A(:,:,i)))*se3_vec(logm(MeanA^(-1)*A(:,:,i)))';
%         SigB = SigB + se3_vec(logm(MeanB^(-1)*B(:,:,i)))*se3_vec(logm(MeanB^(-1)*B(:,:,i)))';
%     end
%     SigA = SigA*(1/size(A,3));
%     SigB = SigB*(1/size(B,3));

    [ VA, ~ ] = eig( SigA(1:3,1:3) );
    [ VB, ~ ] = eig( SigB(1:3,1:3) );

    Q1 = eye(3);
    Q2 = [-1 0 0; 0 -1 0; 0 0 1];
    Q3 = [-1 0 0; 0 1 0; 0 0 -1];
    Q4 = [1 0 0; 0 -1 0; 0 0 -1];
    
    %% 8 solution
    Rx_solved(:,:,1) = VA*Q1*VB';
    Rx_solved(:,:,2) = VA*Q2*VB';
    Rx_solved(:,:,3) = VA*Q3*VB';
    Rx_solved(:,:,4) = VA*Q4*VB';
    Rx_solved(:,:,5) = VA*-Q1*VB';
    Rx_solved(:,:,6) = VA*-Q2*VB';
    Rx_solved(:,:,7) = VA*-Q3*VB';
    Rx_solved(:,:,8) = VA*-Q4*VB';

    %% in fact, just to extract the screw params
    [~, Na, ~, ~] = param_extract(MeanA);
    [~, Nb, ~, ~] = param_extract(MeanB);
    %% rotation axis
    na = so3_vec(Na);
    nb = so3_vec(Nb);
    
    %% traverse the 8 sols and find the minimum
    min = inf;
    for i = 1:8
        if (abs(det(Rx_solved(:,:,i))-1)<0.001) && (norm(na-Rx_solved(1:3,1:3,i)*nb) < min)
            min = norm(na-Rx_solved(1:3,1:3,i)*nb);
            Rx = Rx_solved(:,:,i);
        else
        end

    end
    
    %% equ 20b
    tx_temp = so3_vec(((Rx'*SigA(1:3,1:3)*Rx)^-1*(SigB(1:3,4:6)-Rx'*SigA(1:3,4:6)*Rx))');
    tx = -Rx*tx_temp;

    X = [Rx tx; [0 0 0] 1];

    t_error = (MeanA(1:3,1:3) - eye(3))*tx - Rx*MeanB(1:3,4) + MeanA(1:3,4);
    t_error = norm(t_error);

end