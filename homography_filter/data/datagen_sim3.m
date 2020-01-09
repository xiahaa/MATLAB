function [P,Hreal,Rs,ws,xi,Etareal] = datagen_sim3(N,dt)
    %% plane
    n0 = [0 0 1]';
    d0 = 1;
    %% point number
    n = 100;
    P = [randn(2,n)*2;1*ones(1,n)];
    %% trajectory
    w = 1/30*pi;
    ts = ((1:N)-1)*dt;
    r = 3; 
    h = 4;
    xi = [r*cos(ts*w);r*sin(ts*w);-h*ones(1,N)];
    xid = [-r*sin(ts*w)*w;r*cos(ts*w)*w;zeros(1,N)];
    %% rotation
    R0 = eye(3);
    wr = randn(3,1)*3;
    Rs = zeros(3,3,N);
    Rs(:,:,1) = R0;
    ws = zeros(3,N);
    ws(:,1) = wr;
    Hreal = zeros(3,3,N);
    
    n = Rs(:,:,1)'*n0;
    d = d0-n0'*xi(:,1);
    
    H1 = Rs(:,:,1)+xi(:,1)*n'/d;
    Hreal(:,:,1) = H1./(det(H1)^(1/3));
    Etareal(:,:,1) = Rs(:,:,1)'*xid(:,1)*n'/d;
    for i = 2:N
        wr(1:2) = randn(2,1).*0.05;
        Rs(:,:,i) = Rs(:,:,i-1)*expm(toso3(wr*dt));
        ws(:,i) = wr;
        n = Rs(:,:,i)'*n0;
        d = d0-n0'*xi(:,i);
        H1 = Rs(:,:,i)+xi(:,i)*n'/d;
        Hreal(:,:,i) = H1./(det(H1)^(1/3));
        Etareal(:,:,i) = Rs(:,:,i)'*xid(:,i)*n'/d;
    end
end