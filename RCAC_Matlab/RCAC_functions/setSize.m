%%
% Compute sizes based on System parameters
lx = size(A,1);
lu = size(B,2);
lz = size(C,1);
ly = size(C,1);

% Set sizes
n   = max(size(pole(sys)));
nu  = max(size(Zs));
d   = 0;

for ii = 1:10
    if sum(C*A^(ii-1)*B) ~= 0
        d = ii;
        break
    end
end
