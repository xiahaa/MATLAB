function varargout = Essential_est_five_point(varargin)
%% This file estimate the essential matrix with the help of the matlab symbolic toolbox
    
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
    varargout{1} = Efinal;
end