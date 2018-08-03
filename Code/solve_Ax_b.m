function Setting =  solve_Ax_b(Setting)

A =  Setting.Lu;
b = Setting.b;
numel = Setting.numel;
uidx = Setting.uidx;
labels = Setting.labels;
seeds = Setting.seeds;
active_label = Setting.active_label;

% LLT Decomposition
[L,p] = chol(A,'lower');
LT = conj(L)';
y = L\b;
x = LT\y;

% Holder for saving to external format
xmat = zeros(numel,1);

for i = 1:size(x,1)
    xmat(uidx(i)) = 1 - x(i);
end

for i = 1:length(seeds)
    if(labels(i)==active_label)
        xmat(seeds(i)) = 1;
    else
        xmat(seeds(i)) = 0;
    end
end

Setting.Xu = xmat;

end