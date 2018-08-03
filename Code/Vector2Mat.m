function Mat = Vector2Mat(Vector,R,C)
%% Vector To Mat
Mat = zeros(R,C);
for k = 1:C
    s = (k-1)*R+1 ;
    e = s + R - 1;
    Mat(:,k) = Vector(s:e);
end
end