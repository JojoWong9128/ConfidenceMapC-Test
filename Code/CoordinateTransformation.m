function [imgOut,R,thetaBase,Rmin,d,ThetaMax] = CoordinateTransformation(img)
[row,col] = size(img);
level = graythresh(img);
BW = im2bw(img,level);
FirstNonZeroCoordiante = zeros(col,2);
for i = 1:col
    index = find(BW(:,i)~=0);
    if length(index)~=0
        FirstNonZeroCoordiante(i,1) = i;
        FirstNonZeroCoordiante(i,2) = index(1);
        FirstNonZeroCoordiante(i,3) = index(end);
    else
        FirstNonZeroCoordiante(i,1) = i;
        FirstNonZeroCoordiante(i,2) = row;
        FirstNonZeroCoordiante(i,3) = 0;
    end
end
BaseLine = round(col/2);
LeftData = FirstNonZeroCoordiante(1:BaseLine,:);
LeftAscend = sort(LeftData(:,2),'ascend');
LeftMin = LeftAscend(1);
LeftMinIndx = find(LeftData(:,2)==LeftMin);
LeftMinX = LeftData(LeftMinIndx,1);
LeftMinY = ones(size(LeftMinX)).*LeftMin;
LeftVertex = [max(LeftMinIndx),LeftMin];

RightData = FirstNonZeroCoordiante(BaseLine:end,:);
RightMin = min(RightData(:,2)); 
RightMinIndx = find(RightData(:,2)==RightMin);
RightMinX = RightData(RightMinIndx,1);
RightMinY = ones(size(RightMinX)).*RightMin;
RightVertex = [min(RightMinX),RightMin];

figure(1),plot(FirstNonZeroCoordiante(:,1),FirstNonZeroCoordiante(:,2),'-k')
hold on
plot(LeftMinX,LeftMinY,'*r')
hold on
plot(RightMinX,RightMinY,'*b')

% --------- First Horizontal Line -------- %
[Ah,Bh,Ch] = GeneralLineExpression(LeftVertex,RightVertex);
hold on 
x = 1:col;
y = (-Ah*x -Ch)/Bh;
plot(x,y,'b-')

% -------- Base Line(Polar Axis) ------- %
Cv = Ah*(LeftVertex(2)+RightVertex(2))/2 - Bh*(LeftVertex(1)+RightVertex(1))/2;
Av = Bh;
Bv = -Ah;
hold on
y = 1:row;
x = (Ah*y-Cv)/Bh;
plot(x,y,'r-');
PolarAxisIntensity = img(y,round(x));
LineInt =  PolarAxisIntensity(:,1);
minY = min(find(LineInt~=0));
maxY = max(find(LineInt~=0));
minYX = (Ah*minY-Cv)/Bh;
maxYX = (Ah*maxY-Cv)/Bh;
d = sqrt((minY-maxY)^2+(minYX-maxYX)^2);

%%  Rmin
VertexDeltaX = RightVertex(1) - LeftVertex(1);
LeftVertex2Right = FirstNonZeroCoordiante(LeftVertex(1):RightVertex(1),:) ;
VDeltaY = max(LeftVertex2Right(:,2)) - (LeftVertex(2) + RightVertex(2))/2;
Rmin = (VDeltaY^2 + (VertexDeltaX/2)^2)/(2*VDeltaY);
%% 2*ThetaMax()
% VertexSumX = RightVertex(1) + LeftVertex(1);
% VertexSumY = LeftVertex(2) + RightVertex(2);
% VertexDeltaY = RightVertex(2) - LeftVertex(2);
a = VertexDeltaX^2 + 1;
b = -(2*RightVertex(1)+2*VertexDeltaX*((RightVertex(1)^2-LeftVertex(1)^2+RightVertex(2)^2-LeftVertex(2)^2)/2-RightVertex(2)));
c = RightVertex(1)^2 + ((RightVertex(1)^2-LeftVertex(1)^2+RightVertex(2)^2-LeftVertex(2)^2)/2-RightVertex(2))^2 -round(Rmin)^2;
rootCondition = b^2-4*a*c;
if rootCondition >= 0
    x0 = (-b+sqrt(rootCondition))/(2*a);
    y0 = (RightVertex(1)^2-LeftVertex(1)^2+RightVertex(2)^2-LeftVertex(2)^2)/2 - VertexDeltaX*x0;
    % -------- theta estimation ------------- %
    innerProduct_V1V2 = (LeftVertex(1)-x0)*(RightVertex(1)-x0)+(LeftVertex(2)-y0)*(RightVertex(2)-y0);
    magV1 = sqrt((LeftVertex(1)-x0)^2+(LeftVertex(2)-y0)^2);
    magV2 = sqrt((RightVertex(1)-x0)^2+(RightVertex(2)-y0)^2);
    CosTheta = innerProduct_V1V2/(magV1*magV2);
    ThetaMax = acos(CosTheta)/(pi)*180;
else
    disp('No root in domain R');
end
%%
% % ------------ Right Boundary ------------ %
% RightLine = FirstNonZeroCoordiante(RightVertex(1):end,:);
% [Ar,Br,Cr] = LineFitting(RightLine(1:end-100,:));
% hold on
% Xr = RightLine(:,1);
% Yest = (-Ar*Xr-Cr)/Br;
% plot(Xr,Yest,'-g')
% % --------- Theta_max Calculation -------- %
% CosTheta = (Ah*Br-Bh*Ar)/(sqrt(Ah^2+Bh^2)*sqrt(Ar^2+Br^2));
% ThetaMax = acos(CosTheta)/(pi)*180;
% % ------------ Left Boundary ------------ %
% % LeftLineP1 = FirstNonZeroCoordiante(1:100,:);
% % LeftLineP2 = FirstNonZeroCoordiante(135:LeftVertex(1),:);
% % LeftLine = [LeftLineP1;LeftLineP2];
% % [Al, Bl, Cl] = LineFitting(LeftLine(2:end,:));
% % hold on
% % Xl = LeftLine(:,1);
% % Yest = (-Al*Xl-Cl)/Bl;
% % plot(Xl,Yest,'-g')
% % ------------- Original Point ----------- %
% D = Av*Br-Ar*Bv;
% x0 = (Bv*Cr - Br*Cv)/D;
% y0 = (Ar*Cv - Av*Cr)/D;
% hold on
% plot(x0,y0,'m*')
% Rmin = sqrt((RightVertex(1)-x0)^2+(RightVertex(2)-y0)^2);
% 
% % ------------- Arc Length --------------- %
% % ShortArc = FirstNonZeroCoordiante(LeftVertex(1):RightVertex(1),:);
% % BottomIndex = find(ShortArc(:,2)==max(ShortArc(:,2)));
% % BottomX = 0;
% % if length(BottomIndex)>1
% %     BottomX = median(ShortArc(BottomIndex,1));
% % else
% %     BottomX = ShortArc(BottomIndex);
% % end
%%
ArcLength_min = ThetaMax*pi*Rmin/180;
% VertexRightX = FirstNonZeroCoordiante(end,1);
% VertexRightY = (-VertexRightX*Ar-Cr)/Br;
% d = sqrt((VertexRightY-RightVertex(2))^2+(VertexRightX-RightVertex(1))^2);
ArcLength_max = ThetaMax*pi*(Rmin+d)/180;
Rmax = round(Rmin + d);
% ---- Polar to Cartesian Coordinate ---- %
newImg = zeros(d,ArcLength_min);
RealPart = ones(2*Rmax-1)*diag(1:2*Rmax-1) - Rmax*ones(2*Rmax-1);
ImagPart = 1i*(diag(1:2*Rmax-1)*ones(2*Rmax-1) - Rmax*ones(2*Rmax-1));
Index = RealPart + ImagPart;
thetaBase = mod(-angle(Index)+pi/2+acos(CosTheta)/2,2*pi);
R = abs(Index);
newR = R;
newR(newR<Rmin)=0;
newR(newR>Rmax)=0;
thetaBase(thetaBase>acos(CosTheta))=0;
thetaBase(newR==0)=0;
R(newR==0)=0;
newHeight = max(FirstNonZeroCoordiante(:,3)) - min(LeftVertex(2),RightVertex(2)) + 1;
[TUVindY,TUVindX] = find(thetaBase~=0);
RS = min(TUVindY)- min(LeftVertex(2),RightVertex(2)) ;
CS = Rmax;
thetaBase = thetaBase(RS:RS+newHeight-1,CS-col/2+1:CS+col/2);
R = R(RS:RS+newHeight-1,CS-col/2+1:CS+col/2);
Base = thetaBase;
Base(Base~=0)=1;
R = R.*Base;
imgSeg = Base.*double(img(1:newHeight,:));

for i = 1:size(imgSeg,1)
    for j = 1:size(imgSeg,2)
        Rho = round(R(i,j));
        Theta = thetaBase(i,j)+pi;
        c = round(ArcLength_max/(2*pi)*Theta);
        if Rho~=0 && c~=0
            newImg(Rho,c) = imgSeg(i,j);
        end
    end
end
[indy,indx] = find(newImg~=0);
imgOut = newImg(min(indy)+5:end,min(indx)+5:end);

end

function  [A,B,C] = GeneralLineExpression(Point1,Point2)
% Ax+By+C=0
% syms A B C;
% [A,B,C] = solve(A*Point1(1) + B*Point1(2) + C == 0,A*Point2(1) + B*Point2(2) + C == 0,A,B,C);
A = Point2(2) - Point1(2);
B = Point1(1) - Point2(1);
C = Point2(1)*Point1(2) - Point1(1)*Point2(2);
end

function [A, B, C] = LineFitting(Data)
% Ax + By + C = 0
% Optimized Method : Least Square
Xmean = mean(Data(:,1));
Ymean = mean(Data(:,2));
DxxTemp = (Data(:,1) - Xmean).*(Data(:,1) - Xmean);
Dxx = sum(DxxTemp);
DyyTemp = (Data(:,2) - Ymean).*(Data(:,2) - Ymean);
Dyy = sum(DyyTemp);
DxyTemp = (Data(:,1) - Xmean).*(Data(:,2) - Ymean);
Dxy = sum(DxyTemp);
lambda = ((Dxx+Dyy)-sqrt((Dxx-Dyy)*(Dxx-Dyy)+4*Dxy*Dxy))/2;
den = sqrt(Dxy*Dxy+(lambda-Dxx)*(lambda-Dxx));
A = Dxy/den;
B = (lambda-Dxx)/den;
C = -A*Xmean - B*Ymean;
end