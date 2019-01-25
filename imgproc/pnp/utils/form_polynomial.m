function out = form_polynomial(p, q)
    % 3 * n
    % q should be normalized
    d12 = norm(p(:,1) - p(:,2));
    d13 = norm(p(:,1) - p(:,3));
    d23 = norm(p(:,2) - p(:,3));

    c1 = norm(q(:,1));
    c2 = norm(q(:,2));
    c3 = norm(q(:,3));
    cos12 = q(:,1)'*q(:,2) / c1 / c2;
    cos13 = q(:,1)'*q(:,3) / c1 / c3;
    cos23 = q(:,2)'*q(:,3) / c2 / c3;
    
    poly13 = [1 1 -2*cos13 -d13*d13];
    poly23 = [1 1 -2*cos23 -d23*d23];
    poly12_4order = sylvester_resultant_22(poly13, poly23);
    
    poly12 = [1 1 -2*cos12 -d12*d12];
    poly1_8order = sylvester_resultant_42(poly12_4order, poly12);
    
    %% flip
    out = fliplr(poly1_8order);
end