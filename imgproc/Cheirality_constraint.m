function valid = Cheirality_constraint(E,R,t,q1,q2,K)
    %% triangulation
    valid = zeros(1,size(R,3));
    for i = 1:size(R,3) 
        R1 = R(:,:,i);t1 = t(:,:,i);
        Q = fast_triangulation(q1,q2,E,R1,t1,K);
        if (Q(3)*Q(4) < 0)
            valid(i) = 0;continue;
        else
            c1 = [R1 t1]*Q;
            if (c1(3)*C(4) < 0)
                valid(i) = 0;continue;
            end
        end
        valid(i) = 1;
    end
end