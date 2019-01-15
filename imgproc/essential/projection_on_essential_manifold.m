function Eh = projection_on_essential_manifold(E)
%% basically, it aims to find a matrix embedded on essential manifold which
% minimizes the Frobnius norm $$ norm(Eh-E)_{frob} $$.
    [U,~,V] = svd(E);
    Eh = U*diag([1,1,0])*V';
end