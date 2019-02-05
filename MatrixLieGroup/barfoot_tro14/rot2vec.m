function [ phi ] = rot2vec( C )
% ROT2VEC Compute the matrix log of the rotation matrix C.
%
% From: Timothy D Barfoot and Paul T Furgale,
%       Associating Uncertainty with Three-Dimensional Poses for use in Estimation Problems
%		DOI: 10.1109/TRO.2014.2298059
%       tim.barfoot@utoronto.ca, paul.furgale@mavt.ethz.ch
%
%		Implementation of algorithm 2b in Appendix B, in the paper.
%
% 		For implementation notes on a better version of this function,
% 		please refer to the following discussion:
% 		http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToAngle/
%
% input:
%   C: a 3x3 rotation matrix
%
% output:
%   phi: a 3x1 vector (axis * angle) computed from C
%

rotValidate(C);
phi = [1;0;0];
[v,d]=eig(C);

cond = (trace(C)-1)/2;
if abs(cond+1) > 1e-10
    for i=1:3
       if abs(d(i,i)-1) < 1e-10
          a = v(:,i);
          a = a/sqrt(a'*a);
          phim = acos((trace(C)-1)/2);
          phi = phim*a;

          if abs(trace(vec2rot( phi )'*C)-3) > 1e-14
             phi = -phi;
          end
       end
    end
else
    logA = zeros(3,3);
    for i = 1:2
        logA = logA + (-1)^(i-1)/i.*(C-eye(3))^(i);
    end
    %% project to so3
    logA = 0.5.*(logA-logA');
    phi = [-logA(2,3);logA(1,3);-logA(1,2)];
end

end


