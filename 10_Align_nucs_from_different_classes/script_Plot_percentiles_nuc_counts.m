Figures_folder = 'Figures';
mkdir(Figures_folder)

%% S2_exp1
load('../../data/S2_exp1/Fit_results.100-200.NF40.S2_exp1.mat')

Quantiles_Nuc_counts = quantile(NucCount, [0.05:0.05:0.95]);

x = [1, 2, 5, 15, 40, 60];
C = [1, .4, .4];
dC = ([1, .85, .85] - C) / 8;

figure('Position', [50, 50, 400, 300]);
hold all
for q = [1,5]
    p(q) = patch([x, fliplr(x)],[Quantiles_Nuc_counts(20-q,:), fliplr(Quantiles_Nuc_counts(q,:))], C+(9-q)*dC, 'EdgeColor', 'none');
end
l = plot(x, Quantiles_Nuc_counts(10,:), 'r', 'linewidth', 2);
    
set(gca, 'FontSize', 14)
xlabel('Digestion time (min)')
ylabel('Nucleosome counts')
h = legend(p([1,5]), {'5%-95%', '25%-75%'}, 'location', 'SW');

v = get(h,'title');
set(v,'string','Percentiles', 'FontSize', 12);
set(gcf, 'PaperPositionMode', 'auto')
print(gcf, '-depsc', '-painters', fullfile(Figures_folder, 'Percentiles_nuc_counts.S2_exp1.eps'));


%% S2_exp2
load('../../data/S2_exp2/Fit_results.100-200.NF40.S2_exp2.mat')

Quantiles_Nuc_counts = quantile(NucCount, [0.05:0.05:0.95]);

x = [1, 2, 5, 15, 40, 60];
C = [1, .4, .4];
dC = ([1, .85, .85] - C) / 8;

figure('Position', [50, 50, 400, 300]);
hold all
for q = [1,5]
    p(q) = patch([x, fliplr(x)],[Quantiles_Nuc_counts(20-q,:), fliplr(Quantiles_Nuc_counts(q,:))], C+(9-q)*dC, 'EdgeColor', 'none');
end
l = plot(x, Quantiles_Nuc_counts(10,:), 'r', 'linewidth', 2);
    
set(gca, 'FontSize', 14)
xlabel('Digestion time (min)')
ylabel('Nucleosome counts')
h = legend(p([1,5]), {'5%-95%', '25%-75%'}, 'location', 'SW');

v = get(h,'title');
set(v,'string','Percentiles', 'FontSize', 12);
set(gcf, 'PaperPositionMode', 'auto')
print(gcf, '-depsc', '-painters', fullfile(Figures_folder, 'Percentiles_nuc_counts.S2_exp2.eps'));


%% Kc167_exp1
load('../../data/Kc167_exp1/Fit_results.100-200.NF40.Kc167_exp1.mat')

Quantiles_Nuc_counts = quantile(NucCount, [0.05:0.05:0.95]);

x = [1, 2, 5, 15, 40, 60];
C = [1, .4, .4];
dC = ([1, .85, .85] - C) / 8;

figure('Position', [50, 50, 400, 300]);
hold all
for q = [1,5]
    p(q) = patch([x, fliplr(x)],[Quantiles_Nuc_counts(20-q,:), fliplr(Quantiles_Nuc_counts(q,:))], C+(9-q)*dC, 'EdgeColor', 'none');
end
l = plot(x, Quantiles_Nuc_counts(10,:), 'r', 'linewidth', 2);
    
set(gca, 'FontSize', 14)
xlabel('Digestion time (min)')
ylabel('Nucleosome counts')
h = legend(p([1,5]), {'5%-95%', '25%-75%'}, 'location', 'SW');

v = get(h,'title');
set(v,'string','Percentiles', 'FontSize', 12);
set(gcf, 'PaperPositionMode', 'auto')
print(gcf, '-depsc', '-painters', fullfile(Figures_folder, 'Percentiles_nuc_counts.Kc167_exp1.eps'));

