function Setting = setImage(vector,ImgHeight,ImgWidth)
Setting.alpha = 2;
normalizeValues = 1;
Setting = setMatrix2D(vector,ImgHeight,ImgWidth,normalizeValues,Setting);
% Mat = Vector2Mat(DistPenalty,ImgHeight,ImgWidth);
% Automatic seed generation for confidence maps
seeds = zeros(ImgWidth*2,1);
labels  = zeros(ImgWidth*2,1);
% Upper row seeds and labels
for i = 1:ImgWidth
    seeds(i) = (i-1)*ImgHeight + 1;
    labels(i) = 0;
end
% Lower row seeds and labels
for i = 1:ImgWidth
    seeds(ImgWidth+i) = (i - 1)*ImgHeight + ImgHeight;
    labels(ImgWidth+i) = 1;
end
Setting = setLabeling(seeds,labels,0,2,Setting);
end

function Setting = setMatrix2D(vector,rows,cols,normalizeValues,Setting)

Setting.vector = [];
Setting.rows = rows;
Setting.cols = cols;
Setting.numel = rows * cols;

if(normalizeValues)
    min_mat = min(vector);
    max_mat = max(vector);
    diff = max_mat - min_mat;
    epsilon_diff = 1.0e-17;
    diff = max(diff,epsilon_diff);
    vector = (vector-min_mat)/diff;
end

dist_penalty = zeros(rows,1);

for i = 1:rows
    dist_penalty(i) = 1.0 - (exp(-Setting.alpha*i/rows));
end

% Apply distance weighting to image
for j = 1:cols
    for i = 2:rows
        vector(i+(j-1)*rows) = vector(i+(j-1)*rows)*dist_penalty(i);
    end
end

for i = 1:rows:Setting.numel - rows + 1
    vector(i) = 0;
end
Setting.vector = vector;
end

function Setting = setLabeling(seeds,labels,active_label,num_labels,Setting)
Setting.seeds = seeds;
Setting.labels = labels;
Setting.active_label = active_label;
Setting.num_labels = num_labels;
Setting.isLaplaceAvailable = false;
end