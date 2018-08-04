% -------- Confidence Map C++ Data Test ---------- %
function main()
clc;clear all;close all

addpath('..\Data');
% ---- Data Loading ------ %
% imgROI = csvread('imgROI_T01.csv');
% [row,col] = size(imgROI(:,1:end-1));

[ImName,ImPath] = uigetfile({'*.png;*.jpg;*.bmp;*.tiff'},'File Selector');

SrcInfo = [ImPath,ImName];

if ~ischar(SrcInfo)
    disp('No image file is selected,please check!')
else
    SrcImg = imread(SrcInfo);
    disp('Please select probe type:');
    disp('1 :  flat probe');
    disp('2 :  convex probe');
%     ProbeType = input('Probe Type:');
    ProbeType = 2;
    
    switch ProbeType
        case 1
            img_width = 457;
            img_height = 477;
            inputImg = rgb2gray(SrcImg(62:62+img_height-1,252:252+img_width-1,:));
        case 2
            disp('crop image or not ?');
            disp('1 :  Yes');
            disp('0 :  No');
%             ImageCropping = input('image-cropping:');
            ImageCropping = 1;
            if ImageCropping
                img = SrcImg(76:530,172:786,:);
                img = MarkerEliminating(img);
            else
                img = rgb2gray(SrcImg);  
            end
            tic
%             addpath('..\..\Data');
%             img = csvread('imgSeg.csv');
            [inputImg,RadiusMatrix,ThetaMatrix, Rmin,d,ThetaMax]= CoordinateTransformation(img);
    end
    Vector_imgROI = Mat2Vector(inputImg);
    [row,col] = size(inputImg);
    %%
    Beta = 90.0;
    Gamma = 0.03;
    %%
    Setting = setImage(Vector_imgROI,row,col);
    Setting = computeMap(Beta,Gamma,Setting); 
    ConfidenceMapVector = Setting.Xu;
    ConfidenceMap = Vector2Mat(ConfidenceMapVector,row,col);
    if ( ProbeType == 2 )
        ConfidenceMapFinal = Cartesian2Polar(ConfidenceMap,RadiusMatrix,ThetaMatrix);
    else
        ConfidenceMapFinal = ConfidenceMap;
    end

    %%
    figure(1),imshow(ConfidenceMapFinal,[])
    figure(2),imshow(SrcImg,[])
end
end



