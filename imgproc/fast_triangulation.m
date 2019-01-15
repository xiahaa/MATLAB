function Qs = fast_triangulation(q1,q2,E,R,t,K)
%% triangulation a batch of points.
    if ~isempty(K)
        q1 = K\q1;
        q2 = K\q2;
    end
    P = [R t];
    Qs = zeros(4,size(q1,2));
    for i = 1:size(q1,2)
        a = E'*q2(:,i);
        b = cross(q1(:,i), diag([1,1,0])*a);
        c = cross(q2(:,i), diag([1,1,0])*E*q1(:,i));
        d = cross(a,b);
        C = P'*c;
        Q = [d'*C(4);-(d(1)*C(1)+d(2)*C(2)+d(3)*C(3))];
        Q = Q./Q(4);
        Qs(:,i) = Q;
    end
end