% Load data
load('../../data/S2_exp1/Fit_results.100-200.NF40.S2_exp1.mat')
[N1,edges] = histcounts(O_vector(O_vector <= 0.999), 10, 'Normalization', 'probability');

% Load data
load('../../data/S2_exp2/Fit_results.100-200.NF35.S2_exp2.mat')
[N2,edges] = histcounts(O_vector(O_vector <= 0.999), 10, 'Normalization', 'probability');

%%
figure
hold all
histogram('BinEdges',edges,'BinCounts',100*N1)
histogram('BinEdges',edges,'BinCounts',100*N2)

%%
mean([100*N1; 100*N2])'