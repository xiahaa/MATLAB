syms a b c d real
M = sym('M',[11,11],'real');
for i = 1:size(M,1)
    M(i:end,i) = M(i,i:end);
end
alpha = [1 a^2 a*b a*c a*d b*b b*c b*d c*c c*d d*d]';
da = diff(alpha,a);
db = diff(alpha,b);
dc = diff(alpha,c);
dd = diff(alpha,d);

f1 = alpha'*M*da;
f2 = alpha'*M*db;
f3 = alpha'*M*dc;
f4 = alpha'*M*dd;

[t1,c1]=coeffs(f1,[a b c d]);
[t2,c2]=coeffs(f2,[a b c d]);
[t3,c3]=coeffs(f3,[a b c d]);
[t4,c4]=coeffs(f4,[a b c d]);

syms l real
A = sym('A',[1,11],'real');
f5 = A*da + 2*a*l;
f6 = A*db + 2*b*l;
f7 = A*dc + 2*c*l;
f8 = A*dd + 2*d*l;


[t5,c5]=coeffs(f5,[a b c d l]);
[t6,c6]=coeffs(f6,[a b c d l]);
[t7,c7]=coeffs(f7,[a b c d l]);
[t8,c8]=coeffs(f8,[a b c d l]);