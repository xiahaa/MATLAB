function varargout = hand_eye_calib_sol()
    
    path = mfilename('fullpath');
    i = findstr(path,'\');
    folder = path(1:i(end));
    cd(folder);

    addpath('./solver/');
    addpath('../MatrixLieGroup/barfoot_tro14');
    addpath('../quaternion');
    
    A1 = [-0.989992 -0.141120 0.000000 0; ...
           0.141120 -0.989992 0.000000 0; ...
           0.000000 0.000000 1.000000 0;
           0 0 0 1];
    B1 = [-0.989992 -0.138307 0.028036 -26.9559; ...
           0.138307 -0.911449 0.387470 -96.1332; ...
          -0.028036 0.387470 0.921456 19.4872; ...
           0 0 0 1];
       
    A2 = [0.070737 0.000000 0.997495 -400.000; ...
          0.000000 1.000000 0.000000  0.000000; ...
          -0.997495 0.000000 0.070737 400.000; ...
          0 0 0 1];
    B2 = [0.070737 0.198172 0.977612 -309.543; ...
         -0.198172 0.963323 -0.180936 59.0244; ...
         -0.977612 -0.180936 0.107415 291.177; ...
          0 0 0 1];
    C = A1(1:3,1:3);
    A1(1:3,1:3) = C * inv(sqrtm(C'*C));
    C = A2(1:3,1:3);
    A2(1:3,1:3) = C * inv(sqrtm(C'*C));
    C = B1(1:3,1:3);
    B1(1:3,1:3) = C * inv(sqrtm(C'*C));
    C = B2(1:3,1:3);
    B2(1:3,1:3) = C * inv(sqrtm(C'*C));

      
    T1(1,:,:)= A1;
    T1(2,:,:)= A2;
    T2(1,:,:)= B1;
    T2(2,:,:)= B2;
    
    N = 10;
    dR = euler2rot(rand(1)*pi-pi*0.5,rand(1)*pi-pi*0.5,rand(1)*2*pi-pi)
    dt = rand(3,1)*5-2.5
    dT = [dR dt;[0 0 0 1]];
    dTt = [dR' -dR'*dt;[0 0 0 1]];
    for i = 1:N
        R1 = euler2rot(rand(1)*pi-pi*0.5,rand(1)*pi-pi*0.5,rand(1)*2*pi-pi);
        t1 = rand(3,1)*20-10;
        TSS = [R1 t1;[0 0 0 1]];
        T1(i,:,:) = TSS;
        T2(i,:,:) = dT * TSS * dTt;
    end
%     tic
%     Ts1 = sol_shiu_admad(T1,T2,2);
%     toc
    tic
    Ts2 = sol_tsai_lenz(T1,T2,N);
    toc
    tic
    Ts3 = sol_park_martin(T1,T2,N);
    toc
    tic
    Ts4 = sol_horaud(T1,T2,N);
    toc
    disp(Ts4*dTt)
    Ts6 = sol_dual_quaternion(T1,T2,N);
    disp(Ts6*dTt)
    Ts9 = sol_improved_dual_quaternion(T1,T2,N);
    disp(Ts9*dTt)
    Ts10 = sol_adjoint_transformation_algo(T1,T2,N);
    disp(Ts10*dTt)
    Ts12 = sol_chou(T1,T2,N);
    disp(Ts12*dTt)
    
    tic
    Ts5 = sol_horaud_nlopt(T1,T2,N);
    Ts5*dTt
    toc
    tic
    Ts7 = sol_cvx2(T1,T2,N);
    disp(Ts7*dTt)
    toc
    tic
    Ts8 = sol_dphec(T1,T2,N);
    disp(Ts8*dTt)
    toc
    tic
    Ts11 = sol_dual_sdp_cvx(T1,T2,N);
    disp(Ts11*dTt)
    toc
    tic
    Ts12 = sol_manifold_opt_SE3(T1,T2,N);
    toc
    disp(Ts12*dTt)
end

