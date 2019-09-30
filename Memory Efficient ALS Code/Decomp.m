function Z = Decomp(Size, j, i)
%DECOMP Compute the kronecker decomposition for a given j and i
Z = [zeros((j-1)*Size,Size); speye(Size,Size); zeros((i+1-j)*Size,Size)];
end

