function img = MarkerEliminating(img)
ChannelR = img(:,:,1);
ChannelG = img(:,:,2);
ChannelB = img(:,:,3);
% figure(1),imshow(ChannelR);
% figure(2),imshow(ChannelG);
% figure(3),imshow(ChannelB);
Marker = ChannelB - ChannelR;
GrayImg = rgb2gray(img).*(uint8(1-logical(Marker)));
targetArea1 = GrayImg(300:400,580:610);
maxInt = max(max(GrayImg));
[indX1,indY1] = find(targetArea1 == maxInt);
if length(indX1)>4
    h = fspecial('average',11);
    filteredImg = imfilter(targetArea1,h);
    targetArea1(indX1,indY1) = filteredImg(indX1,indY1);

    targetArea2 = GrayImg(:,570:end);
    [indX2,indY2] = find(targetArea2 >= maxInt*0.9);
    targetArea2(indX2,indY2)=0;
    GrayImg(:,570:end) = targetArea2;
    GrayImg(300:400,580:610) = targetArea1;
end
img = GrayImg;
end