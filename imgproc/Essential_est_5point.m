function E = Essential_est_5point(varargin)
%% Essential matrix estimation using exactly 5 point.
% see David Niester: An efficient solution to the five-point relative pose
% estimation problem. IEEE PAMI.
    p1 = varargin{1};
    p2 = varargin{2};
    K = varargin{3};
    Kinv = inv(K);
    p1n = Kinv*p1;
    p2n = Kinv*p2;
    
    A = zeros(5,9);
    A(1,:) = [p1n(:,1)'*p2n(1,1) p1n(:,1)'*p2n(2,1) p1(:,1)'*p2n(3,1)];
    A(2,:) = [p1n(:,2)'*p2n(1,2) p1n(:,2)'*p2n(2,2) p1(:,2)'*p2n(3,2)];
    A(3,:) = [p1n(:,3)'*p2n(1,3) p1n(:,3)'*p2n(2,3) p1(:,3)'*p2n(3,3)];
    A(4,:) = [p1n(:,4)'*p2n(1,4) p1n(:,4)'*p2n(2,4) p1(:,4)'*p2n(3,4)];
    A(5,:) = [p1n(:,5)'*p2n(1,5) p1n(:,5)'*p2n(2,5) p1(:,5)'*p2n(3,5)];
    [~,~,V] = svd(A);
    mask = [1 2 3;4 5 6;7 8 9];
    X = V(:,end-3);
    Y = V(:,end-2);
    Z = V(:,end-1);
    W = V(:,end);
    X = X(mask);Y = Y(mask);Z = Z(mask);W = W(mask);
    
    d11 = constraint1(X,Y,Z,W);
    d12 = constaint2(X,Y,Z,W);
    [Agj,~,~] = formA(d11,d12);
    xyz = xyzest(Agj);
    
    E = zeros(3,3,size(xyz,2));
    for i = 1:1:size(xyz,2)
        x = xyz(1,i);y = xyz(2,i);z = xyz(3,i);
        E(:,:,i) = x*X+y*Y+z*Z+W;
        E(:,:,i) = E(:,:,i) ./ E(3,3,i);
%         diag(p2'*Eest(:,:,i)*p1)
    end
end

function varargout = xyzest(Agj)
% xzz xz x yzz yz y zzz zz z 1
%  1  2   3   4   5  6   7  8   9 10 11 12 13
    e = Agj(5,11:20);
    f = Agj(6,11:20);
    g = Agj(7,11:20);
    h = Agj(8,11:20);
    i = Agj(9,11:20);
    j = Agj(10,11:20);
    
    k = [0 e(1:3) 0 e(4:6) 0 e(7:10)] - [f(1:3) 0 f(4:6) 0 f(7:10) 0];
    l = [0 g(1:3) 0 g(4:6) 0 g(7:10)] - [h(1:3) 0 h(4:6) 0 h(7:10) 0];
    m = [0 i(1:3) 0 i(4:6) 0 i(7:10)] - [j(1:3) 0 j(4:6) 0 j(7:10) 0];
    
    B = cell(3,3);
    B{1,1} = k(1:4);B{1,2} = k(5:8);B{1,3} = k(9:13);
    B{2,1} = l(1:4);B{2,2} = l(5:8);B{2,3} = l(9:13);
    B{3,1} = m(1:4);B{3,2} = m(5:8);B{3,3} = m(9:13);
    
    %% method 1
%     tic
%     detB = z334(B{1,1},B{2,2},B{3,3}) + z334(B{1,2},B{3,1},B{2,3}) + z334(B{3,2},B{2,1},B{1,3}) ...
%          - z334(B{3,1},B{2,2},B{1,3}) - z334(B{1,2},B{2,1},B{3,3}) - z334(B{1,1},B{3,2},B{2,3});
% %     toc
%     zs = roots(detB); 
%     N = numel(find(abs(imag(zs)) < 1e-8));
%     xyz = zeros(3,N);
%     k = 1;
%     for i = 1:10
%         if abs(imag(zs(i))) < 1e-8
%             %% real root
%             z = real(zs(i));
%             A = [peval(B{1,1},z),peval(B{1,2},z),peval(B{1,3},z); ...
%                  peval(B{2,1},z),peval(B{2,2},z),peval(B{2,3},z); ...
%                  peval(B{3,1},z),peval(B{3,2},z),peval(B{3,3},z)];
%             [~,~,V] = svd(A);
%             sol = V(:,end);
%             xyz(:,k) = [sol(1)/sol(3);sol(2)/sol(3);z];
%             k = k + 1;
%         end
%     end
%     toc
    %% method 2
%     tic
    p1 = z3z4(B{1,2},B{2,3}) - z3z4(B{1,3},B{2,2});
    p2 = z3z4(B{1,3},B{2,1}) - z3z4(B{1,1},B{2,3});
    p3 = z3z3(B{1,1},B{2,2}) - z3z3(B{1,2},B{2,1});
    detB = z7z3(p1,B{3,1}) + z7z3(p2,B{3,2}) + z6z4(p3,B{3,3});
    zs = roots(detB); 
    N = numel(find(abs(imag(zs)) < 1e-8));
    xyz = zeros(3,N);
    k = 1;
    for i = 1:10
        if abs(imag(zs(i))) < 1e-8
            %% real root
            z = real(zs(i));
            den = peval(p3,z);
            x = peval(p1,z)/den;
            y = peval(p2,z)/den;
            xyz(:,k) = [x;y;z];
            k = k + 1;
        end
    end
%     toc
    
    varargout{1} = xyz;
end

function res = z6z4(c1,c2)
% t = [ z^10, z^9, z^8, z^7, z^6, z^5, z^4, z^3, z^2, z, 1]
    if numel(c1) == 5 && numel(c2) == 7
        c11 = c1(1);c12 = c1(2);c13 = c1(3);c14 = c1(4);c15 = c1(5);
        c21 = c2(1);c22 = c2(2);c23 = c2(3);c24 = c2(4);c25 = c2(5);c26 = c2(6);c27 = c2(7);
    elseif numel(c1) == 7 && numel(c2) == 5
        c11 = c2(1);c12 = c2(2);c13 = c2(3);c14 = c2(4);c15 = c2(5);
        c21 = c1(1);c22 = c1(2);c23 = c1(3);c24 = c1(4);c25 = c1(5);c26 = c1(6);c27 = c1(7);
    end
    res = [ c11*c21, ...
            c11*c22 + c12*c21, ...
            c11*c23 + c12*c22 + c13*c21, ...
            c11*c24 + c12*c23 + c13*c22 + c14*c21, ...
            c11*c25 + c12*c24 + c13*c23 + c14*c22 + c15*c21, ...
            c11*c26 + c12*c25 + c13*c24 + c14*c23 + c15*c22, ...
            c11*c27 + c12*c26 + c13*c25 + c14*c24 + c15*c23, ...
            c12*c27 + c13*c26 + c14*c25 + c15*c24, ...
            c13*c27 + c14*c26 + c15*c25, ...
            c14*c27 + c15*c26, ...
            c15*c27];
end

function res = z7z3(c1,c2)
%% t = [ z^10, z^9, z^8, z^7, z^6, z^5, z^4, z^3, z^2, z, 1]
    if numel(c1) == 4 && numel(c2) == 8
        c11 = c1(1);c12 = c1(2);c13 = c1(3);c14 = c1(4);
        c41 = c2(1);c42 = c2(2);c43 = c2(3);c44 = c2(4);c45 = c2(5);c46 = c2(6);c47 = c2(7);c48 = c2(8);
    elseif numel(c1) == 8 && numel(c2) == 4
        c11 = c2(1);c12 = c2(2);c13 = c2(3);c14 = c2(4);
        c41 = c1(1);c42 = c1(2);c43 = c1(3);c44 = c1(4);c45 = c1(5);c46 = c1(6);c47 = c1(7);c48 = c1(8);
    end
    res = [ c11*c41, ...
            c11*c42 + c12*c41, ...
            c11*c43 + c12*c42 + c13*c41, ...
            c11*c44 + c12*c43 + c13*c42 + c14*c41, ...
            c11*c45 + c12*c44 + c13*c43 + c14*c42, ...
            c11*c46 + c12*c45 + c13*c44 + c14*c43, ...
            c11*c47 + c12*c46 + c13*c45 + c14*c44, ...
            c11*c48 + c12*c47 + c13*c46 + c14*c45, ...
            c12*c48 + c13*c47 + c14*c46, ...
            c13*c48 + c14*c47, ...
            c14*c48];

end

function res = z3z3(c1,c2)
% t = [ z^6, z^5, z^4, z^3, z^2, z, 1]
    c11 = c1(1);c12 = c1(2);c13 = c1(3);c14 = c1(4);
    c21 = c2(1);c22 = c2(2);c23 = c2(3);c24 = c2(4);
    res = [ c11*c21, ...
            c11*c22 + c12*c21, ...
            c11*c23 + c12*c22 + c13*c21, ...
            c11*c24 + c12*c23 + c13*c22 + c14*c21, ...
            c12*c24 + c13*c23 + c14*c22, ...
            c13*c24 + c14*c23, ...
            c14*c24];
end

function res = z3z4(c1,c2)
%% [ z^7, z^6, z^5, z^4, z^3, z^2, z, 1]
    if numel(c1) == 4 && numel(c2) == 5
        c11 = c1(1);c12 = c1(2);c13 = c1(3);c14 = c1(4);
        c31 = c2(1);c32 = c2(2);c33 = c2(3);c34 = c2(4);c35 = c2(5);
    elseif numel(c1) == 5 && numel(c2) == 4
        c11 = c2(1);c12 = c2(2);c13 = c2(3);c14 = c2(4);
        c31 = c1(1);c32 = c1(2);c33 = c1(3);c34 = c1(4);c35 = c1(5);
    end

    res = [ c11*c31, ...
            c11*c32 + c12*c31, ...
            c11*c33 + c12*c32 + c13*c31, ...
            c11*c34 + c12*c33 + c13*c32 + c14*c31, ...
            c11*c35 + c12*c34 + c13*c33 + c14*c32, ...
            c12*c35 + c13*c34 + c14*c33, ...
            c13*c35 + c14*c34, ...
            c14*c35];
end

function res = peval(c, x)
    res = 0;
    n = numel(c);
    for i = 1:numel(c)
        res = res + c(i) * x^(n-i);
    end
end

function res = z334(c1,c2,c3)
%% [ z^10, z^9, z^8, z^7, z^6, z^5, z^4, z^3, z^2, z, 1]        
    c11 = c1(1);c12 = c1(2);c13 = c1(3);c14 = c1(4);
    c21 = c2(1);c22 = c2(2);c23 = c2(3);c24 = c2(4);
    c31 = c3(1);c32 = c3(2);c33 = c3(3);c34 = c3(4);c35 = c3(5);
    
    res = [ c11*c21*c31, ...
            c31*(c11*c22 + c12*c21) + c11*c21*c32, ...
            c32*(c11*c22 + c12*c21) + c31*(c11*c23 + c12*c22 + c13*c21) + c11*c21*c33, ...
            c33*(c11*c22 + c12*c21) + c32*(c11*c23 + c12*c22 + c13*c21) + c31*(c11*c24 + c12*c23 + c13*c22 + c14*c21) + c11*c21*c34, ...
            c34*(c11*c22 + c12*c21) + c33*(c11*c23 + c12*c22 + c13*c21) + c31*(c12*c24 + c13*c23 + c14*c22) + c32*(c11*c24 + c12*c23 + c13*c22 + c14*c21) + c11*c21*c35, ...
            c35*(c11*c22 + c12*c21) + c31*(c13*c24 + c14*c23) + c34*(c11*c23 + c12*c22 + c13*c21) + c32*(c12*c24 + c13*c23 + c14*c22) + c33*(c11*c24 + c12*c23 + c13*c22 + c14*c21), ...
            c32*(c13*c24 + c14*c23) + c35*(c11*c23 + c12*c22 + c13*c21) + c33*(c12*c24 + c13*c23 + c14*c22) + c34*(c11*c24 + c12*c23 + c13*c22 + c14*c21) + c14*c24*c31, c33*(c13*c24 + c14*c23) + c34*(c12*c24 + c13*c23 + c14*c22) + c35*(c11*c24 + c12*c23 + c13*c22 + c14*c21) + c14*c24*c32, ...
            c34*(c13*c24 + c14*c23) + c35*(c12*c24 + c13*c23 + c14*c22) + c14*c24*c33, ...
            c35*(c13*c24 + c14*c23) + c14*c24*c34, ...
            c14*c24*c35];
end

function varargout = formA(d1,d2)
% [ x^3, x^2*y, x^2*z, x^2, x*y^2, x*y*z, x*y, x*z^2, x*z, x, y^3, y^2*z, 
%     1      2      3    4      5      6    7      8    9  10  11     12
%   y^2, y*z^2, y*z, y, z^3, z^2, z, 1]
%    13     14   15 16   17   18  19 20
    Afull = [d1;d2];% 10x20
    mask = [1 11 2 5 3 4 12 13 6 7 8 9 10 14 15 16 17 18 19 20];
    A = Afull(:,mask);
    
    [~,U] = lu(A);%% LU decomposition
    Agj = U;
    for i = size(Agj,1):-1:5
        Agj(i,:) = Agj(i,:)./Agj(i,i);%% pivoting value to 1
        for j = i+1:1:size(A,1)
            Agj(i,:) = Agj(i,:) - Agj(j,:).*Agj(i,j);
        end
    end
    varargout{1} = Agj;
    varargout{2} = A;
    varargout{3} = Afull;
end

function d1 = constraint1(X,Y,Z,W)
    p11 = [X(1,1),Y(1,1),Z(1,1),W(1,1)];
    p12 = [X(1,2),Y(1,2),Z(1,2),W(1,2)];
    p13 = [X(1,3),Y(1,3),Z(1,3),W(1,3)];
    
    p21 = [X(2,1),Y(2,1),Z(2,1),W(2,1)];
    p22 = [X(2,2),Y(2,2),Z(2,2),W(2,2)];
    p23 = [X(2,3),Y(2,3),Z(2,3),W(2,3)];
    
    p31 = [X(3,1),Y(3,1),Z(3,1),W(3,1)];
    p32 = [X(3,2),Y(3,2),Z(3,2),W(3,2)];
    p33 = [X(3,3),Y(3,3),Z(3,3),W(3,3)];
    
    d1 = p2p1(p1p1(p12,p23)-p1p1(p13,p22),p31) + ...
         p2p1(p1p1(p13,p21)-p1p1(p11,p23),p32) + ...
         p2p1(p1p1(p11,p22)-p1p1(p12,p21),p33);
end

function d2 = constaint2(X,Y,Z,W)
    p11 = [X(1,1),Y(1,1),Z(1,1),W(1,1)];
    p12 = [X(1,2),Y(1,2),Z(1,2),W(1,2)];
    p13 = [X(1,3),Y(1,3),Z(1,3),W(1,3)];
    
    p21 = [X(2,1),Y(2,1),Z(2,1),W(2,1)];
    p22 = [X(2,2),Y(2,2),Z(2,2),W(2,2)];
    p23 = [X(2,3),Y(2,3),Z(2,3),W(2,3)];
    
    p31 = [X(3,1),Y(3,1),Z(3,1),W(3,1)];
    p32 = [X(3,2),Y(3,2),Z(3,2),W(3,2)];
    p33 = [X(3,3),Y(3,3),Z(3,3),W(3,3)];


    EEt11 = p1p1(p11,p11) + p1p1(p12,p12) + p1p1(p13,p13);
    EEt12 = p1p1(p11,p21) + p1p1(p12,p22) + p1p1(p13,p23);
    EEt13 = p1p1(p11,p31) + p1p1(p12,p32) + p1p1(p13,p33);
    EEt22 = p1p1(p21,p21) + p1p1(p22,p22) + p1p1(p23,p23);
    EEt23 = p1p1(p21,p31) + p1p1(p22,p32) + p1p1(p23,p33);
    EEt33 = p1p1(p31,p31) + p1p1(p32,p32) + p1p1(p33,p33);
    
    C1 = 0.5.*(EEt11+EEt22+EEt33);
    Delta11 = EEt11 - C1;
    Delta22 = EEt22 - C1;
    Delta33 = EEt33 - C1;
    
    Delta12 = EEt12;Delta13 = EEt13;Delta23 = EEt23;
    Delta21 = Delta12;Delta31 = Delta13;Delta32 = Delta23;
    
    d2 = zeros(9,20);
    d2(1,:) = p2p1(Delta11,p11) + p2p1(Delta12,p21) + p2p1(Delta13,p31);% 11
    d2(2,:) = p2p1(Delta11,p12) + p2p1(Delta12,p22) + p2p1(Delta13,p32);
    d2(3,:) = p2p1(Delta11,p13) + p2p1(Delta12,p23) + p2p1(Delta13,p33);
    
    d2(4,:) = p2p1(Delta21,p11) + p2p1(Delta22,p21) + p2p1(Delta23,p31);
    d2(5,:) = p2p1(Delta21,p12) + p2p1(Delta22,p22) + p2p1(Delta23,p32);
    d2(6,:) = p2p1(Delta21,p13) + p2p1(Delta22,p23) + p2p1(Delta23,p33);
    
    d2(7,:) = p2p1(Delta31,p11) + p2p1(Delta32,p21) + p2p1(Delta33,p31);
    d2(8,:) = p2p1(Delta31,p12) + p2p1(Delta32,p22) + p2p1(Delta33,p32);
    d2(9,:) = p2p1(Delta31,p13) + p2p1(Delta32,p23) + p2p1(Delta33,p33);
end

function res = p1p1(c1,c2)
%% assume order like x y z 1, so polyc1 * polyc2 = 
% [ x^2, x*y, x*z, x, y^2, y*z, y, z^2, z, 1]
    c11 = c1(1);c12 = c1(2);c13 = c1(3);c14 = c1(4);
    c21 = c2(1);c22 = c2(2);c23 = c2(3);c24 = c2(4);

    res = [ c11*c21, ...
            c11*c22 + c12*c21, ...
            c11*c23 + c13*c21, ...
            c11*c24 + c14*c21, ...
            c12*c22, ...
            c12*c23 + c13*c22, ...
            c12*c24 + c14*c22, ...
            c13*c23, ...
            c13*c24 + c14*c23, ...
            c14*c24];
end

function res = p2p1(c1,c2)
% [ x^3, x^2*y, x^2*z, x^2, x*y^2, x*y*z, x*y, x*z^2, x*z, x, y^3, y^2*z, y^2, y*z^2, y*z, y, z^3, z^2, z, 1]
    p11 = c1(1);p12 = c1(2);p13 = c1(3);p14 = c1(4);p15 = c1(5);p16 = c1(6);p17 = c1(7);p18 = c1(8);p19 = c1(9);p10 = c1(10);
    c11 = c2(1);c12 = c2(2);c13 = c2(3);c14 = c2(4);
    
    res = [ c11*p11, ...
            c11*p12 + c12*p11, ...
            c11*p13 + c13*p11, ...
            c11*p14 + c14*p11, ...
            c12*p12 + c11*p15, ...
            c12*p13 + c13*p12 + c11*p16, ...
            c12*p14 + c14*p12 + c11*p17, ...
            c13*p13 + c11*p18, ...
            c13*p14 + c14*p13 + c11*p19, ...
            c11*p10 + c14*p14, ...
            c12*p15, ...
            c12*p16 + c13*p15, ...
            c12*p17 + c14*p15, ...
            c13*p16 + c12*p18, ...
            c13*p17 + c14*p16 + c12*p19, ...
            c12*p10 + c14*p17, ...
            c13*p18, ...
            c13*p19 + c14*p18, ...
            c13*p10 + c14*p19, ...
            c14*p10];
end