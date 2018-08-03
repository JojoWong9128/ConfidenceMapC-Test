function Vector = Mat2Vector(Mat)
[row,col] = size(Mat);
Vector = zeros(row*col,1);
for k = 1:col
    s = (k-1)*row+1;
    e = s + row - 1;
    Vector(s:e) = Mat(:,k);
end
end