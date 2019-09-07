function xi = data_term_error(Rdata,Rreg,indices,varargin)
    if nargin == 3
        xi = zeros(3,length(indices));
        for i = 1:length(indices)
            ii = indices(i);
            xi(:,i) = logSO3(Rdata(:,(i*3-2):i*3)'*Rreg(:,(ii*3-2):ii*3));
            xi(:,i) = para_trans(Rdata(:,(i*3-2):i*3),Rreg(:,(ii*3-2):ii*3),xi(:,i));
        end
    else
        id = find(indices == varargin{1},1);
        if isempty(id) 
            xi = [];
        else
            ii = indices(id);
            xi = logSO3(Rdata(:,(id*3-2):id*3)'*Rreg(:,(ii*3-2):ii*3));
            xi = para_trans(Rdata(:,(id*3-2):id*3),Rreg(:,(ii*3-2):ii*3),xi);
        end
    end
end