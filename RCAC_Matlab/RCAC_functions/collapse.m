function [Acol] = collapse(A)

    m       = size(A,1);
    n       = size(A,2);
    k       = size(A,3);
    Acol    = zeros(m*k,n);

    for ii = 1:k
        Acol(m*(ii-1)+1:m*ii,:) = A(:,:,ii);
    end
end