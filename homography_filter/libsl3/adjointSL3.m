function AdH = adjointSL3(H,X)
    AdH = H*X*inv(H);
end