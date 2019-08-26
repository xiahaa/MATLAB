
function vpt = para_trans(R1,R2,v)
    A = logSO3(R1'*R2);
    Rpt = expSO3(A/2);
    vpt = invhat(R2'*R1*Rpt*hat(v)*Rpt);
end