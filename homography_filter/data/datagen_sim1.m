function [Hreal, Hmeas, Areal] = datagen_sim1(N,dt)
    %% fist simulation
    H0 = [1 0 0;0 1 0;0 0 1];
    A = [0.1 0.1 0.1;0.05 -0.18 0.025;0.05 0 0.08];
    Hreal = zeros(3,3,N);
    Hmeas = zeros(3,3,N);
    Areal = zeros(3,3,N);

    Hreal(:,:,1) = H0 ./ (det(H0)^(1/3));
    Hmeas(:,:,1) = H0 ./ (det(H0)^(1/3));
    Areal(:,:,1) = A;
    for i = 2:N
        Q = randn(3,3).*1;
        Q = tosl3(Q);
        H = Hreal(:,:,i-1)*expm((A).*dt);% exp(Ad_H U dt)*H = H*exp(U*dt), so this equals to expm(Hreal(:,:,i-1)*(A).*dt*inv(Hreal(:,:,i-1)))*Hreal(:,:,i-1)
        Hreal(:,:,i) = H ./ (det(H)^(1/3));
        H = Hreal(:,:,i-1)*expm((A+Q).*dt);% exp(Ad_H U dt)*H = H*exp(U*dt), so this equals to expm(Hreal(:,:,i-1)*(A).*dt*inv(Hreal(:,:,i-1)))*Hreal(:,:,i-1)
        Hmeas(:,:,i) = H ./ (det(H)^(1/3));
        Areal(:,:,i) = A;
    end
end