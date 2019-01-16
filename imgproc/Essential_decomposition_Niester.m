function [R,t] = Essential_decomposition_Niester(E)
    [U,~,V] = svd(E);
    D = [0 1 0;-1 0 0;0 0 1];
    
    R = zeros(3,3,4);
    t = zeros(3,1,4);
    
    if (det(U)<0)
        U(:,3) = -U(:,3);
    end
    if (det(V)<0)
        V(:,3) = -V(:,3);
    end
    
    R(:,:,1) = U*D*V';t(:,:,1) = U(:,3);
    R(:,:,2) = R(:,:,1);t(:,:,2) = -U(:,3);
    R(:,:,3) = U*D'*V';t(:,:,3) = U(:,3);
    R(:,:,4) = R(:,:,3);t(:,:,4) = -U(:,3);
end