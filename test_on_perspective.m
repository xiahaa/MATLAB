function test_on_perspective

R = eye(3,3);
t = [0;0;0];
A = [1000 0 200;0 1000 300;0 0 1];
P = A*[R t];

ptsrc = rand(3,1000).*10 - 5;
ptsrchomo = [ptsrc;ones(1,size(ptsrc,2))];

pt2dhomo = P*ptsrchomo;

pt2d(1,:) = pt2dhomo(1,:)./pt2dhomo(3,:);
pt2d(2,:) = pt2dhomo(2,:)./pt2dhomo(3,:);

AA = zeros(3*size(pt2d,2),12);

for i = 1:1:size(pt2d,2)
    X = ptsrc(1,i);   Y = ptsrc(2,i);  Z = ptsrc(3,i);
    u = pt2d(1,i);   v = pt2d(2,i); 
   
    id = (i - 1)*3;
    
    AA(id+1,1)  = 0; AA(id+1,2)  = -X; AA(id+1,3)  = v*X;
    AA(id+1,4)  = 0; AA(id+1,5)  = -Y; AA(id+1,6)  = v*Y;
    AA(id+1,7) = 0;  AA(id+1,8)  = -Z; AA(id+1,9)  = v*Z;
    AA(id+1,10)  = 0; AA(id+1,11) = -1; AA(id+1,12) = v;

    AA(id+2,1) = X; AA(id+2,2) =  0; AA(id+2,3) = -u*X;
    AA(id+2,4) = Y; AA(id+2,5) =  0; AA(id+2,6) = -u*Y;
    AA(id+2,7) = Z; AA(id+2,8) =  0; AA(id+2,9) = -u*Z;
    AA(id+2,10) = 1; AA(id+2,11) = 0; AA(id+2,12) = -u;

    AA(id+3,1) = -v*X; AA(id+3,2) = u*X; AA(id+3,3) = 0;
    AA(id+3,4) = -v*Y; AA(id+3,5) = u*Y; AA(id+3,6) = 0;
    AA(id+3,7) = -v*Z; AA(id+3,8) = u*Z; AA(id+3,9) = 0;
    AA(id+3,10) = -v;   AA(id+3,11) = u; AA(id+3,12) = 0;
end

AAT = AA'*AA;

[U,S,V] = svd(AAT);

return

