function newImg = Cartesian2Polar(img,RadiusMatrix,ThetaMatrix)
h = size(img,1);
s=size(img);
newImg = zeros(size(RadiusMatrix));
thetaMax = max(ThetaMatrix(:));
uniqueR = unique(RadiusMatrix(:));
Rmin = uniqueR(2);
% index=1i*(diag(1:(2*h-1))*ones((2*h-1))-h*ones(2*h-1))+ones(2*h-1)*diag(1:(2*h-1))-h*ones(2*h-1); % 表示成了复数矩阵, 方便求极角(辐角).
% theta=mod(-angle(index)-pi/2,pi/4); % 求极角(辐角).
% r=abs(index); % 求极径, 以图片b的中心为坐标原点.
for k=1:size(RadiusMatrix,1)
    for m=1:size(RadiusMatrix,2)
        if RadiusMatrix(k,m)-Rmin <= h && RadiusMatrix(k,m) ~= 0
            newImg(k,m,:)=img(min(h,floor(RadiusMatrix(k,m)-Rmin)+1),max(floor(s(2)*ThetaMatrix(k,m)/thetaMax),1),:); % 新的图像b的每一点和原图像a的关系, 是由原图像的点得到的.
        else
            newImg(k,m,:)= 0;  
        end
    end
end

end