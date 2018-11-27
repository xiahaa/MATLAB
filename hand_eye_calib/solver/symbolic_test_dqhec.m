function symbolic_test_dqhec
%     syms q1 q2 q3 q4 t1 t2 t3 real
%     q = [q1;q2;q3;q4];
%     R = q2r(q);
%     TX = [R [t1;t2;t3];[0 0 0 1]];
%     
%     RA = sym('RA', [3,3], 'real');
%     RB = sym('RB', [3,3], 'real');
%     TA = sym('TA', [3,1], 'real');
%     TB = sym('TB', [3,1], 'real');
%     TA = [RA TA;[0 0 0 1]];
%     TB = [RB TB;[0 0 0 1]];
%     
%     M = TA*TX-TX*TB;
%     poly = trace(M'*M);
%     [tt,cc] = coeffs(poly,[q1 q2 q3 q4 t1 t2 t3]);


    %% this is how the dqhec generate, so there is a trick
%     syms qx1 qx2 qx3 qx4 qx5 qx6 qx7 qx8 real
%     syms a1 a2 a3 a4 a5 a6 a7 a8 real
%     syms b1 b2 b3 b4 b5 b6 b7 b8 real
%     
%     dqx = [qx1 qx2 qx3 qx4 qx5 qx6 qx7 qx8]';
%     dqa = [a1 a2 a3 a4 a5 a6 a7 a8]';
%     dqb = [b1 b2 b3 b4 b5 b6 b7 b8]';
%     l2 = proddualquaternion(dqx,dqb);
%     dqxc = conjugatedualquaternion(dqx);
%     l3 = proddualquaternion(l2,dqxc);     
%     ld = dqa - l3;    
%     poly = ((ld'*ld));
    

    %% new formulation
    syms q1 q2 q3 q4 t1 t2 t3 real
    q = [q1;q2;q3;q4];
    qt = [0;t1;t2;t3];
    qd = 0.5.*q2m_left(qt)*q;
    dqx = [q;qd];
    
    syms a1 a2 a3 a4 a5 a6 a7 a8 real
    syms b1 b2 b3 b4 b5 b6 b7 b8 real
    dqa = [a1 a2 a3 a4 a5 a6 a7 a8]';
    dqb = [b1 b2 b3 b4 b5 b6 b7 b8]';
    l1 = proddualquaternion(dqa,dqx);
    l2 = proddualquaternion(dqx,dqb);
    ld = l1 - l2;
    poly = ((ld'*ld));
    
    [tt,cc] = coeffs(poly,[q1 q2 q3 q4 t1 t2 t3]);
    
end

function ddp = proddualquaternion(dq1,dq2)
    ddp14 = q2m_left(dq1(1:4))*dq2(1:4);
    ddp58 = q2m_left(dq1(1:4))*dq2(5:8) + q2m_left(dq1(5:8))*dq2(1:4);
    ddp = [ddp14;ddp58];
end

function ddqc = conjugatedualquaternion(dq)
    ddqc1 = [dq(1);-dq(2:4)];
    ddqc2 = [dq(5);-dq(6:8)];
    ddqc = [ddqc1;ddqc2];
end