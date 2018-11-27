function  varargout = sol_dual_sdp_cvx(TA,TB,N)
%% interface for calling the hand eye calibration methods implemented by 
% Matthew Giamou at University of Toronto
% Certifiably Globally Optimal Extrinsic Calibration from Per-Sensor Egomotion
% ICRA 2018

    dim = size(TA,2);
    if N < 2
        error('At least two samples needed for unique solution!');
        varargout{1} = [];
        return;
    end
    if dim < 4
        error('Only work for T!');
        varargout{1} = [];
        return;
    end
    
    format short;

    R1 = zeros(3,3,N);
    R2 = zeros(3,3,N);
    t1 = zeros(3,N);
    t2 = zeros(3,N);
    
    for i = 1:N
        T1 = TA(i,:,:);T1 = reshape(T1,dim,dim,1);
        T2 = TB(i,:,:);T2 = reshape(T2,dim,dim,1);
        
        R1(:,:,i) = T1(1:3,1:3);
        t1(:,i) = T1(1:3,4);
        R2(:,:,i) = T2(1:3,1:3);
        t2(:,i) = T2(1:3,4);        
    end

    addpath('C:/Users/xiahaa/Documents/DTU/paper/hand_eye_calibration/code/certifiable-calibration-master/certifiable-calibration-master/matlab/egomotion_calibration');
    addpath('C:/Users/xiahaa/Documents/DTU/paper/hand_eye_calibration/code/certifiable-calibration-master/certifiable-calibration-master/matlab/utils');

    disp('Estimate:');
    [R_cal, t_cal] = egomotion_calibration(R1,t1,R2,t2)
    disp('Ground Truth:');
    
    if dim == 4
        varargout{1} = [R_cal t_cal;[0 0 0 1]];
    else
        varargout{1} = R12;
    end
end