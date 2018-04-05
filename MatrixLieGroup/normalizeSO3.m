function Rn = normalizeSO3(R)
Rn = R * inv(sqrtm(R'*R));
end

