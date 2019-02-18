
function relative_SO3s = form_relative_SO3(SO3s)
    nSO3s = size(SO3s,3);
    set_r_SO3s = nchoosek(1:nSO3s, 2);
    n_r_SO3s = size(set_r_SO3s,1);
    relative_SO3s = zeros(3,3,n_r_SO3s);
    for i = 1:1:n_r_SO3s
        ii = set_r_SO3s(i,1);
        jj = set_r_SO3s(i,2);
        relative_SO3s(:,:,i) = SO3s(:,:,ii)*SO3s(:,:,jj)';
    end
end