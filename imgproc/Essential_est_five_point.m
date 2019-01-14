function varargout = Essential_est_five_point(varargin)

%     clear;
%     clc;
%     randn('seed',0);
%     R = [0.4120    0.1411    0.9002;
%         0.7969   -0.5349   -0.2808;
%         0.4419    0.8330   -0.3328];
%     tran = [1 2 3]';
%     Etruth = [0 -3 2;
%               3 0 -1;
%               -2 1 0]*R;
%     Etruth = Etruth/Etruth(3,3)
% 
%     p1a = randn(3,1);
%     p1b = R*p1a+tran;
%     p1a = p1a/p1a(3);
%     p1b = p1b/p1b(3);
% 
%     p2a = randn(3,1);
%     p2b = R*p2a+tran;
%     p2a = p2a/p2a(3);
%     p2b = p2b/p2b(3);
% 
%     p3a = randn(3,1);
%     p3b = R*p3a+tran;
%     p3a = p3a/p3a(3);
%     p3b = p3b/p3b(3);
% 
%     p4a = randn(3,1);
%     p4b = R*p4a+tran;
%     p4a = p4a/p4a(3);
%     p4b = p4b/p4b(3);
% 
%     p5a = randn(3,1);
%     p5b = R*p5a+tran;
%     p5a = p5a/p5a(3);
%     p5b = p5b/p5b(3);
% 
%     p6a = randn(3,1);
%     p6b = R*p6a+tran;
%     p6a = p6a/p6a(3);
%     p6b = p6b/p6b(3);
% 
%     p1 = [p1a p2a p3a p4a p5a];
%     p2 = [p1b p2b p3b p4b p5b];
    
    p1 = varargin{1};
    p2 = varargin{2};
    
    if size(p1,2) ~= 5
        error('5 points !!');
    end
    
    %% 1st step
    p1 = p1 ./ p1(3,:);
    p2 = p2 ./ p2(3,:);
    A = zeros(5,9);
    for i = 1:5
        A(i,:) = [p1(:,i)'.*p2(1,i) p1(:,i)'.*p2(2,i) p1(:,i)'.*p2(3,i)];
    end
    [~,~,V] = svd(A);
    X = V(:,end-3);
    Y = V(:,end-2);
    Z = V(:,end-1);
    W = V(:,end);
    
    mask = [1 2 3;4 5 6;7 8 9];
    X = X(mask);Y = Y(mask);Z = Z(mask);W = W(mask);
    
    %% symbolic
    syms x y z real
    E = x.*X + y.*Y + z.*Z + W;
    d1 = (det(E));
    d2 = (E*E'*E - 0.5.*trace(E*E')*E);
    d = [d1;d2(:)];
        
    Cz = sym('coeffd',[10,10],'real');
    for i = 1:10
        [Cz(i,:), ~] = coeffs(d(i), [x y]);
    end
    detcz = det(Cz);
    [coeffdcz, ~] = coeffs(detcz, z);
    
    %% roots
    rootsz = roots(coeffdcz);
    
    validid = (abs(imag(eval(rootsz))) < 1e-8);
    xyz = zeros(3,numel(find(validid)));
    k = 1;
    for i = 1:size(rootsz,1)
        if validid(i) == 1
            Czz = subs(Cz, rootsz(i));
            [~,~,V] = svd(Czz);
            xyz(:,k) = [V(6,end)/V(10,end);V(9,end)/V(10,end);rootsz(i)];
            k = k + 1;
        end
    end
    Efinal = zeros(3,3,size(xyz,2));
    for i = 1:size(xyz,2)
        Efinal(:,:,i) = xyz(1,i).*X + xyz(2,i).*Y + xyz(3,i).*Z + W;
        Efinal(:,:,i) = Efinal(:,:,i) ./ Efinal(3,3,i); 
%         diag(p2'*Efinal(:,:,i)*p1)
    end
end