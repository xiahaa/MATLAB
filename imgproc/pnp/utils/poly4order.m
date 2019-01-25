function varargout = poly4order(p,qn)
    cos12 = qn(:,1)'*qn(:,2);
    cos13 = qn(:,1)'*qn(:,3);
    cos23 = qn(:,2)'*qn(:,3);
    
    l0 = 1;l1 = cos12;l2 = cos13;
    C1 = sin(acos(cos12));
    C2 = sin(acos(cos13));
    
    D1 = norm(p(:,1)-p(:,2));
    D2 = norm(p(:,1)-p(:,3));
    D3 = norm(p(:,2)-p(:,3));
    
    if abs(D1) < 1e-6
        varargout{1} = [];
        return;
    end
    
    k = D2/D1;
    A1 = k*k;
    A2 = A1*C1*C1-C2*C2;
    A3 = l2*cos23-l1;
    A4 = l1*cos23-l2;
    A5 = cos23;
    D12 = D1*D1;
    A6 = (D3*D3-D12-D2*D2)/(2*D12);
    A7 = 1 - l1*l1 - l2*l2 + l1*l2*cos23 + A6*C1*C1;
    
    B4 = A6*A6-A1*A5*A5;
    B3 = 2*(A3*A6-A1*A4*A5);
    B2 = A3*A3+2*A6*A7-A1*A4*A4-A2*A5*A5;
    B1 = 2*(A3*A7-A2*A4*A5);
    B0 = A7*A7-A2*A4*A4;
    
    varargout{1} = [B4 B3 B2 B1 B0];
    varargout{2} = [A1 A2 A3 A4 A5 A6 A7];
    varargout{3} = D1;
    varargout{4} = [l1 l2];
    varargout{5} = C1;
end