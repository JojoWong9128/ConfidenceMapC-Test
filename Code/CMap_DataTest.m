% -------- Confidence Map C++ Data Test ---------- %
function main()
clc;clear all;close all

addpath('..\Data');
% ---- Data Loading ------ %
imgROI = csvread('imgROI_T01.csv');
[row,col] = size(imgROI(:,1:end-1));
% [row,col] = size(imgROI);
Vector_imgROI = Mat2Vector(imgROI);
%%
Beta = 90.0;
Gamma = 0.03;
%%
Setting = setImage(Vector_imgROI,row,col);
Setting = computeMap(Beta,Gamma,Setting); 
MapMatrix = Vector2Mat(Setting.Xu,row,col);
%%
figure(1),imshow(MapMatrix,[])
figure(2),imshow(imgROI,[])
end



