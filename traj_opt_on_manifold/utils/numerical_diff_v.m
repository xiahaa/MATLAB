function v = numerical_diff_v(Rreg,varargin)
    N = round(size(Rreg,2)/3);
    if nargin == 1
        v = zeros(3,N-1);
        for i = 1:N-1
            v(:,i) = logSO3(Rreg(:,(i*3-2):i*3)'*Rreg(:,(i*3+1):(i*3+3)));
        end
    else
        i = varargin{1};
        if i == 1
            % tangent plane 2
            v(:,1) = logSO3(Rreg(:,(i*3+1):i*3+3)'*Rreg(:,(i*3-2):(i*3)));% 2 -- 1            
            v(:,2) = logSO3(Rreg(:,(i*3+1):i*3+3)'*Rreg(:,(i*3+4):(i*3+6)));% 2 -- 3
            
            v(:,1) = para_trans(Rreg(:,(i*3+1):i*3+3),Rreg(:,(i*3-2):i*3),v(:,1));% tgt2 to 1
            v(:,2) = para_trans(Rreg(:,(i*3+1):i*3+3),Rreg(:,(i*3-2):i*3),v(:,2));% tgt2 to 1

            % old code
%             v(:,1) = logSO3(Rreg(:,(i*3-2):i*3)'*Rreg(:,(i*3+1):(i*3+3)));
%             v(:,2) = logSO3(Rreg(:,(i*3+1):i*3+3)'*Rreg(:,(i*3+4):(i*3+6)));
        elseif i == N
            v(:,1) = logSO3(Rreg(:,(i*3-5):i*3-3)'*Rreg(:,(i*3-8):(i*3-6)));%e-1 -- e-2
            v(:,2) = logSO3(Rreg(:,(i*3-5):i*3-3)'*Rreg(:,(i*3-2):(i*3)));%e-1 -- e
            
            v(:,1) = para_trans(Rreg(:,(i*3-5):i*3-3),Rreg(:,(i*3-2):i*3),v(:,1));% e-1 -- e
            v(:,2) = para_trans(Rreg(:,(i*3-5):i*3-3),Rreg(:,(i*3-2):i*3),v(:,2));
            
%             v(:,1) = logSO3(Rreg(:,(i*3-8):i*3-6)'*Rreg(:,(i*3-5):(i*3-3)));
%             % tangent plane n-1
%             v(:,2) = logSO3(Rreg(:,(i*3-5):i*3-3)'*Rreg(:,(i*3-2):(i*3)));
%             
        elseif i == 2
            % little trick, 1--2 parallel to 2 equals to -(1--2)
            v(:,1) = -logSO3(Rreg(:,(i*3-2):(i*3))'*Rreg(:,(i*3-5):(i*3-3)));% 1 -- 2 pt to 2
            v(:,2) = -logSO3(Rreg(:,(i*3-2):(i*3))'*Rreg(:,(i*3+1):(i*3+3)));% 3 -- 2 pt to 2
            v(:,3) = logSO3(Rreg(:,(i*3+1):i*3+3)'*Rreg(:,(i*3+4):(i*3+6)));% 3 -- 4
            
            v(:,3) = para_trans(Rreg(:,(i*3+1):i*3+3),Rreg(:,(i*3-2):i*3),v(:,3));% pt to 2
        elseif i == N-1
            v(:,1) = logSO3(Rreg(:,(i*3-5):i*3-3)'*Rreg(:,(i*3-8):(i*3-6)));% e-2 -- e-3
            v(:,2) = -logSO3(Rreg(:,(i*3-2):i*3)'*Rreg(:,(i*3-5):(i*3-3)));% e-2 -- e-1 to e-1
            v(:,3) = -logSO3(Rreg(:,(i*3-2):i*3)'*Rreg(:,(i*3+1):(i*3+3)));% e -- e-1 to e-1
            
            v(:,1) = para_trans(Rreg(:,(i*3-5):(i*3-3)),Rreg(:,(i*3-2):i*3),v(:,1));% pt from e-2 to e-1
        else
            v(:,1) = logSO3(Rreg(:,(i*3-5):i*3-3)'*Rreg(:,(i*3-8):(i*3-6)));% i-1 -- i-2
            v(:,2) = logSO3(Rreg(:,(i*3-5):i*3-3)'*Rreg(:,(i*3-2):(i*3))); % i-1 -- i
            v(:,3) = logSO3(Rreg(:,(i*3+1):i*3+3)'*Rreg(:,(i*3-2):(i*3))); % i+1 -- i
            v(:,4) = logSO3(Rreg(:,(i*3+1):i*3+3)'*Rreg(:,(i*3+4):(i*3+6))); % i+1 -- i+2
            
            v(:,1) = para_trans(Rreg(:,(i*3-5):(i*3-3)),Rreg(:,(i*3-2):i*3),v(:,1));% pt from i-1 to i
            v(:,2) = para_trans(Rreg(:,(i*3-5):(i*3-3)),Rreg(:,(i*3-2):i*3),v(:,2));% pt from i-1 to i
            v(:,3) = para_trans(Rreg(:,(i*3+1):(i*3+3)),Rreg(:,(i*3-2):i*3),v(:,3));% pt from i+1 to i
            v(:,4) = para_trans(Rreg(:,(i*3+1):(i*3+3)),Rreg(:,(i*3-2):i*3),v(:,4));% pt from i+1 to i
        end
    end
end