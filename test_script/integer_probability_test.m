function integer_probability_test

q = 1;
n = -40:1:40;
[S1,n] = getIntegerProbability(n,q)
q = 2;
[S2,n] = getIntegerProbability(n,q)
q=4;
[S3,n] = getIntegerProbability(n,q)


plot(n,S1,'b+');
hold on;
plot(n,S2,'rx');
plot(n,S3,'ko');


return

function [S,n] = getIntegerProbability(n,q)
if (size(n,2) ~= 1)
    n = n';
end
S = zeros(size(n,1),1);
negid = n<=0;
posid = n>0;
S(negid) = 1./(2-q.*n(negid));
S(posid) = (q.*n(posid)+1)./(2+q.*n(posid));
return