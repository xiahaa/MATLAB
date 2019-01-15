function E = Essential_est_8_point(q1, q2)
    %% it is assume that q1 and q2 are 3xN array. N is assumed to be 8.
    A = [q1'.*q2(1,:)' q1'.*q2(2,:)' q1'.*q2(3,:)'];
    [~,~,V] = svd(A);
    E = V(:,end);
    mask = [1 2 3;4 5 6;7 8 9];
    E = E(mask);
end