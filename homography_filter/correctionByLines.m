function w = correctionByLines(Hest, lref, lcur)
% be very careful here, Hest transforms current feature to refered
% feature. Another concern is the |e|.
%
    e = inv(Hest)'*lcur;
    e = e ./ sqrt(e(1,:).^2+e(2,:).^2+e(3,:).^2);
    w = zeros(3,3);
    for i = 1:size(e,2)
        w = w + e(:,i)*lref(:,i)'*(eye(3) - e(:,i)*e(:,i)');
    end
end
