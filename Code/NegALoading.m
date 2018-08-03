function matNegA = NegALoading(NegA)
[row,col] = size(NegA);
% maxRow = max(NegA(:,1)) + 1;
% maxCol = max(NegA(:,2)) + 1;
matNegA = sparse(477,457);
for i = 1:row
    idx = NegA(i,1) + 1;
    idy = NegA(i,2) + 1;
    matNegA(idx,idy) = NegA(i,3);
end
end