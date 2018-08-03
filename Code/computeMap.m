function Setting = computeMap(beta,gamma,Setting)
addpath('..\Data');
mapVector = zeros(Setting.rows*Setting.cols,1);
Setting.beta = beta;
Setting.gamma = gamma;
Setting = solve(Setting);
mapVector = Setting.Xu;
end

function Setting = solve(Setting)

if(isempty(Setting.vector) || isempty(Setting.seeds) || isempty(Setting.labels))
    disp('ERROR: External inputs (vector, seeds, labels, etc.) not available.');
else
    if(~Setting.isLaplaceAvailable)
        tic
        Setting = assembleLaplacianEdges(Setting);  % -A
        % NegA = csvread('NegA.csv');
        % Setting.L = NegALoading(NegA);
        
%         SparseL = SparseMatrixLoading();
%         Setting.L = SparseL;
        toc
        tic
        Setting = assembleDegree(Setting);          % D
        toc
        Setting = generateUniqueIndices(Setting);
        Setting.isLaplaceAvailable = true;
    end
    Setting = assemble_Lu_b(Setting);
    Setting =  solve_Ax_b(Setting);
end

end

function Setting = assembleLaplacianEdges(Setting)

% Numerical limit to avoid division by zero and add to the weights to avoid zero weights in Laplacian
epsilon = 1.0e-5;
L = sparse(Setting.numel,Setting.numel);

min_weight = Setting.vector(1);
max_weight = Setting.vector(1);

% Horizontal edges
for i = 1:Setting.numel - Setting.rows 
    weight = abs(Setting.vector(i) - Setting.vector(i+Setting.rows));
    weight = weight + Setting.gamma;     
    
    if weight < min_weight
        min_weight = weight;
    end
    if weight > max_weight
        max_weight = weight;
    end
end
% Vertical edges
for i = 1:Setting.rows: Setting.numel - Setting.rows + 1
    for j = 1:Setting.rows-1
        weight = abs(Setting.vector(i+j-1) - Setting.vector(i+j));
        if weight < min_weight
            min_weight = weight;
        end
        if weight > max_weight
            max_weight = weight;
        end
    end
end

diff = (max_weight - min_weight);
epsilon_diff = 1.0e-17;

if diff < epsilon_diff
    diff = epsilon_diff;
end
% L.reserve(VectorXi::Constant(numel,5)); 缺省这一步
% VectorXi::Constant(numel,5)即：5*ones(numel,1);

for i = 1:Setting.numel - Setting.rows 
    weight = abs(Setting.vector(i) - Setting.vector(i+Setting.rows));
    weight = (weight - min_weight)/diff;
    weight = weight + Setting.gamma; 
    weight = exp(-Setting.beta*weight)+epsilon;
    
    L(i,i+Setting.rows) = -weight;
    L(i+Setting.rows,i) = -weight;
end
for i = 1:Setting.rows: Setting.numel - Setting.rows + 1
    for j = 1:Setting.rows-1
        weight = abs(Setting.vector(i+j -1) - Setting.vector(i+j));
        weight = (weight - min_weight)/diff;
        weight = exp(-Setting.beta*weight)+epsilon;
        L(i+j-1,i+j) = -weight;
        L(i+j,i+j-1) = -weight;
    end    
end
Setting.L = L;
end

function Setting = assembleDegree(Setting)
for i = 1:Setting.numel
    Setting.L(i,i) = abs(sum(Setting.L(:,i)));
end
end

function Setting = generateUniqueIndices(Setting)
n = length(Setting.seeds);   % marked nodes
q = Setting.numel - n;       % unmarked nodes

% uidx_temp = zeros(Setting.numel,1);
% Set all node indices
% for i = 1:Setting.numel
%     uidx_temp(i) = i;
% end

 uidx_temp = linspace(1,Setting.numel,Setting.numel); 
 
 % Flag marked node with negative sign (last -1 incase of 0th node being a seed)
 for i = 1:n
     uidx_temp(Setting.seeds(i)) = -uidx_temp(Setting.seeds(i))-1;
 end
 % Vector of unmarked nodes. Marked nodes are in seeds vector
 uidx = zeros(q,1);
 k = 0;
 for i = 1:Setting.numel
     if(uidx_temp(i)>=0)
         k = k+1;
         uidx(k) = uidx_temp(i);
     end
 end
Setting.uidx = uidx;
end

function Setting = assemble_Lu_b(Setting)
n = length(Setting.seeds);   % marked nodes
q = Setting.numel - n;       % unmarked nodes
uidx = Setting.uidx;
Cu = sparse(Setting.numel,q);

for i = 1:q
    Cu(uidx(i),i) = 1;
end
% Cu.transpose() means actually rows unmarked. Once again A^T C A is encountered.
Lu = Cu'*Setting.L*Cu;
% Permutation matrix to get MARKED COLUMNS
Cm = sparse(Setting.numel,n);

for i = 1:n
    Cm(Setting.seeds(i),i) = 1;
end

Bt = sparse(q,n);
Bt = Cu'*Setting.L*Cm;

M = sparse(n,Setting.num_labels);

for i = 1 : Setting.num_labels
    for j = 1 : n
        if(Setting.labels(j) == i) 
            M(j,i) = 1;
        end
    end
end

CL = sparse(size(M,2),1);
CL(1) = 1;
Setting.b = -Bt*M*CL;
Setting.Lu = Lu;
end