function Q = fast_triangulation(q1,q2,E,P)
    a = E'*q2;
    b = cross(q1, diag([1,1,0])*a);
    c = cross(q2, diag([1,1,0])*E*q1);
    d = cross(a,b);
    C = P'*c;
    Q = [d'*C(4);-(d(1)*C(1)+d(2)*C(2)+d(3)*C(3))];
    Q = Q./Q(4);
end