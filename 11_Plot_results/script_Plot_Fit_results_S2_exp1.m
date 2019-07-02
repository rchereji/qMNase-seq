% Load data
load('../../data/S2_exp1/Fit_results.100-200.NF40.S2_exp1.mat')
norm_factor = 40;

Figures_folder = 'S2_exp1';
mkdir(Figures_folder)

%%
figure('Position',[100 100 400 300]);
Counts = histc(k1_vector, -0.01:0.02:50.01);
Counts = Counts(1:end-1);
BinCenters = 0:0.02:50;

bar(BinCenters, Counts, 1)
xlim([0 50])
ylim([0 21000])
xlabel('Release rate constant (k1)')
ylabel('Number of nucleosomes')
set(gca, 'FontSize', 14)
set(gcf, 'PaperPositionMode', 'auto')
print(gcf, '-depsc', '-painters', fullfile(Figures_folder, 'k1_dist_linear_scale.eps'));

P = get(gca, 'Position');

%%
set(gca, 'Position', [0.2    0.15    0.5    0.5]);
set(gca,'yscale','log')
ylabel({'Number of nucleosomes';'(log scale)'})
xlim([0 50])
print(gcf, '-depsc', '-painters', fullfile(Figures_folder, 'k1_dist_log_scale.eps'));


%%
nuc_id = find(k1_vector >= 10);
myChr = Chr(nuc_id);
myLoc = Loc(nuc_id);
save(fullfile(Figures_folder, 'Fragile_complexes_Dmel.mat'), 'myChr', 'myLoc')

chrReference = {'2L'; '2R'; '3L'; '3R'; '4'; 'X'; 'Y'; 'rDNA'};
chrLabel = chrReference(myChr);

BEDcontents = [];

noNucs = numel(myChr);
for n = 1:noNucs
    BEDcontents = [BEDcontents, sprintf('%s\t%d\t%d\tfc_%d\t0\t*\t%d\t%d\t100,100,100\n', ...
                             chrLabel{n}, myLoc(n)-73, myLoc(n)+73, n, myLoc(n)-73, myLoc(n)+73)];                         
end

% Save BED file
fileID = fopen(fullfile(Figures_folder, 'Dmel_fragile_complexes.bed'), 'w');
fprintf(fileID, '%s', BEDcontents);
fclose(fileID);

%%
nuc_id = find(k1_vector >= 10);

ratio = norm_factor / max(sum(NucCount, 2));
RescalledNucCount = ratio * NucCount;

Quantiles_Nuc_counts = quantile(RescalledNucCount(nuc_id, :), [0.05:0.05:0.95]);
noQuantiles = size(Quantiles_Nuc_counts, 1);

x = [1, 2, 5, 15, 40, 60];
C = [1, .4, .4];
dC = ([1, .85, .85] - C) / 8;

f = figure('Position', [50, 50, 450, 250]);
hold all
for q = [1,5]
    p(q) = patch([x, fliplr(x)],[Quantiles_Nuc_counts(20-q,:), fliplr(Quantiles_Nuc_counts(q,:))], C+(9-q)*dC, 'EdgeColor', 'none');
end
l = plot(x, Quantiles_Nuc_counts(10,:), 'r', 'linewidth', 2);
    
set(gca, 'FontSize', 11)
xlabel('Digestion time (min)', 'FontSize', 12)
ylabel('Apparent occupancy', 'FontSize', 12)
title('Data', 'FontSize', 12)
h = legend(p([1,5]), {'5%-95%', '25%-75%'}, 'location', 'EO');

v = get(h,'title');
set(v,'string','Percentiles', 'FontSize', 11);

xlim([0 15])
ylim([0 1])
set(gca, 'FontSize', 14)
set(gcf, 'PaperPositionMode', 'auto')
print(gcf, '-depsc', '-painters', fullfile(Figures_folder, 'Percentiles_nuc_counts_high_k1.eps'));

% Plot fitted curves
no_sites = length(nuc_id);
F = @(x, m) x(1) * x(2)/(x(2)-(x(3)))*(exp(-(x(3))*m) - exp(-x(2)*m));

FittedOcc = nan(no_sites, 1001);
m = 0:0.1:100;
for s = 1:no_sites
    FittedOcc(s,:) = F([O_vector(nuc_id(s)), k1_vector(nuc_id(s)), k2_vector(nuc_id(s))], m);
end

Quantiles_Nuc_counts = quantile(FittedOcc, [0.05:0.05:0.95]);
noQuantiles = size(Quantiles_Nuc_counts, 1);

x = m;
C = [1, .4, .4];
dC = ([1, .85, .85] - C) / 8;

f = figure('Position', [50, 50, 450, 250]);
hold all
for q = [1,5]
    p(q) = patch([x, fliplr(x)],[Quantiles_Nuc_counts(20-q,:), fliplr(Quantiles_Nuc_counts(q,:))], C+(9-q)*dC, 'EdgeColor', 'none');
end
l = plot(x, Quantiles_Nuc_counts(10,:), 'r', 'linewidth', 2);

set(gca, 'FontSize', 11)
xlabel('Digestion time (min)', 'FontSize', 12)
ylabel('Apparent occupancy', 'FontSize', 12)
title('Fitted', 'FontSize', 12)
h = legend(p([1,5]), {'5%-95%', '25%-75%'}, 'location', 'EO');

v = get(h,'title');
set(v,'string','Percentiles', 'FontSize', 11);

xlim([0 15])
ylim([0 1])
set(gca, 'FontSize', 14)
set(gcf, 'PaperPositionMode', 'auto')
print(gcf, '-depsc', '-painters', fullfile(Figures_folder, 'Percentiles_nuc_counts_high_k1_fitted.eps'));


%%
nuc_id = find((k1_vector < 10) & (O_vector < 0.2));

ratio = norm_factor / max(sum(NucCount, 2));
RescalledNucCount = ratio * NucCount;

Quantiles_Nuc_counts = quantile(RescalledNucCount(nuc_id, :), [0.05:0.05:0.95]);
noQuantiles = size(Quantiles_Nuc_counts, 1);

x = [1, 2, 5, 15, 40, 60];
C = [1, .4, .4];
dC = ([1, .85, .85] - C) / 8;

f = figure('Position', [50, 50, 450, 250]);
hold all
for q = [1,5]
    p(q) = patch([x, fliplr(x)],[Quantiles_Nuc_counts(20-q,:), fliplr(Quantiles_Nuc_counts(q,:))], C+(9-q)*dC, 'EdgeColor', 'none');
end
l = plot(x, Quantiles_Nuc_counts(10,:), 'r', 'linewidth', 2);
    
set(gca, 'FontSize', 11)
xlabel('Digestion time (min)', 'FontSize', 12)
ylabel('Apparent occupancy', 'FontSize', 12)
title('Data; O < 0.2', 'FontSize', 12)
h = legend(p([1,5]), {'5%-95%', '25%-75%'}, 'location', 'EO');

v = get(h,'title');
set(v,'string','Percentiles', 'FontSize', 11);

xlim([0 60])
ylim([0 1])
set(gca, 'FontSize', 14)
set(gcf, 'PaperPositionMode', 'auto')
print(gcf, '-depsc', '-painters', fullfile(Figures_folder, 'Percentiles_nuc_counts_low_k1_bin_1.eps'));

% Plot fitted curves
no_sites = length(nuc_id);
F = @(x, m) x(1) * x(2)/(x(2)-(x(3)))*(exp(-(x(3))*m) - exp(-x(2)*m));

FittedOcc = nan(no_sites, 1001);
m = 0:0.1:100;
for s = 1:no_sites
    FittedOcc(s,:) = F([O_vector(nuc_id(s)), k1_vector(nuc_id(s)), k2_vector(nuc_id(s))], m);
end

Quantiles_Nuc_counts = quantile(FittedOcc, [0.05:0.05:0.95]);
noQuantiles = size(Quantiles_Nuc_counts, 1);

x = m;
C = [1, .4, .4];
dC = ([1, .85, .85] - C) / 8;

f = figure('Position', [50, 50, 450, 250]);
hold all
for q = [1,5]
    p(q) = patch([x, fliplr(x)],[Quantiles_Nuc_counts(20-q,:), fliplr(Quantiles_Nuc_counts(q,:))], C+(9-q)*dC, 'EdgeColor', 'none');
end
l = plot(x, Quantiles_Nuc_counts(10,:), 'r', 'linewidth', 2);

set(gca, 'FontSize', 11)
xlabel('Digestion time (min)', 'FontSize', 12)
ylabel('Apparent occupancy', 'FontSize', 12)
title('Fitted; O < 0.2', 'FontSize', 12)
h = legend(p([1,5]), {'5%-95%', '25%-75%'}, 'location', 'EO');

v = get(h,'title');
set(v,'string','Percentiles', 'FontSize', 11);

xlim([0 60])
ylim([0 1])
set(gca, 'FontSize', 14)
set(gcf, 'PaperPositionMode', 'auto')
print(gcf, '-depsc', '-painters', fullfile(Figures_folder, 'Percentiles_nuc_counts_low_k1_bin_1_fitted.eps'));


%%
nuc_id = find((k1_vector < 10) & (O_vector >= 0.2) & (O_vector < 0.3));

ratio = norm_factor / max(sum(NucCount, 2));
RescalledNucCount = ratio * NucCount;

Quantiles_Nuc_counts = quantile(RescalledNucCount(nuc_id, :), [0.05:0.05:0.95]);
noQuantiles = size(Quantiles_Nuc_counts, 1);

x = [1, 2, 5, 15, 40, 60];
C = [1, .4, .4];
dC = ([1, .85, .85] - C) / 8;

f = figure('Position', [50, 50, 450, 250]);
hold all
% for q = 1:9
for q = [1,5]
    p(q) = patch([x, fliplr(x)],[Quantiles_Nuc_counts(20-q,:), fliplr(Quantiles_Nuc_counts(q,:))], C+(9-q)*dC, 'EdgeColor', 'none');
end
l = plot(x, Quantiles_Nuc_counts(10,:), 'r', 'linewidth', 2);
    
set(gca, 'FontSize', 11)
xlabel('Digestion time (min)', 'FontSize', 12)
ylabel('Apparent occupancy', 'FontSize', 12)
title('Data; 0.2 < O < 0.3', 'FontSize', 12)
h = legend(p([1,5]), {'5%-95%', '25%-75%'}, 'location', 'EO');

v = get(h,'title');
set(v,'string','Percentiles', 'FontSize', 11);

xlim([0 60])
ylim([0 1])
set(gca, 'FontSize', 14)
set(gcf, 'PaperPositionMode', 'auto')
print(gcf, '-depsc', '-painters', fullfile(Figures_folder, 'Percentiles_nuc_counts_low_k1_bin_2.eps'));

% Plot fitted curves
no_sites = length(nuc_id);
F = @(x, m) x(1) * x(2)/(x(2)-(x(3)))*(exp(-(x(3))*m) - exp(-x(2)*m));

FittedOcc = nan(no_sites, 1001);
m = 0:0.1:100;
for s = 1:no_sites
    FittedOcc(s,:) = F([O_vector(nuc_id(s)), k1_vector(nuc_id(s)), k2_vector(nuc_id(s))], m);
end

Quantiles_Nuc_counts = quantile(FittedOcc, [0.05:0.05:0.95]);
noQuantiles = size(Quantiles_Nuc_counts, 1);

x = m;
C = [1, .4, .4];
dC = ([1, .85, .85] - C) / 8;

f = figure('Position', [50, 50, 450, 250]);
hold all
for q = [1,5]
    p(q) = patch([x, fliplr(x)],[Quantiles_Nuc_counts(20-q,:), fliplr(Quantiles_Nuc_counts(q,:))], C+(9-q)*dC, 'EdgeColor', 'none');
end
l = plot(x, Quantiles_Nuc_counts(10,:), 'r', 'linewidth', 2);

set(gca, 'FontSize', 11)
xlabel('Digestion time (min)', 'FontSize', 12)
ylabel('Apparent occupancy', 'FontSize', 12)
title('Fitted; 0.2 < O < 0.3', 'FontSize', 12)
h = legend(p([1,5]), {'5%-95%', '25%-75%'}, 'location', 'EO');

v = get(h,'title');
set(v,'string','Percentiles', 'FontSize', 11);

xlim([0 60])
ylim([0 1])
set(gca, 'FontSize', 14)
set(gcf, 'PaperPositionMode', 'auto')
print(gcf, '-depsc', '-painters', fullfile(Figures_folder, 'Percentiles_nuc_counts_low_k1_bin_2_fitted.eps'));


%%
nuc_id = find((k1_vector < 10) & (O_vector >= 0.3) & (O_vector < 0.4));

ratio = norm_factor / max(sum(NucCount, 2));
RescalledNucCount = ratio * NucCount;

Quantiles_Nuc_counts = quantile(RescalledNucCount(nuc_id, :), [0.05:0.05:0.95]);
noQuantiles = size(Quantiles_Nuc_counts, 1);

x = [1, 2, 5, 15, 40, 60];
C = [1, .4, .4];
dC = ([1, .85, .85] - C) / 8;

f = figure('Position', [50, 50, 450, 250]);
hold all
for q = [1,5]
    p(q) = patch([x, fliplr(x)],[Quantiles_Nuc_counts(20-q,:), fliplr(Quantiles_Nuc_counts(q,:))], C+(9-q)*dC, 'EdgeColor', 'none');
end
l = plot(x, Quantiles_Nuc_counts(10,:), 'r', 'linewidth', 2);
    
set(gca, 'FontSize', 11)
xlabel('Digestion time (min)', 'FontSize', 12)
ylabel('Apparent occupancy', 'FontSize', 12)
title('Data; 0.3 < O < 0.4', 'FontSize', 12)
h = legend(p([1,5]), {'5%-95%', '25%-75%'}, 'location', 'EO');

v = get(h,'title');
set(v,'string','Percentiles', 'FontSize', 11);

xlim([0 60])
ylim([0 1])
set(gca, 'FontSize', 14)
set(gcf, 'PaperPositionMode', 'auto')
print(gcf, '-depsc', '-painters', fullfile(Figures_folder, 'Percentiles_nuc_counts_low_k1_bin_3.eps'));

% Plot fitted curves
no_sites = length(nuc_id);
F = @(x, m) x(1) * x(2)/(x(2)-(x(3)))*(exp(-(x(3))*m) - exp(-x(2)*m));

FittedOcc = nan(no_sites, 1001);
m = 0:0.1:100;
for s = 1:no_sites
    FittedOcc(s,:) = F([O_vector(nuc_id(s)), k1_vector(nuc_id(s)), k2_vector(nuc_id(s))], m);
end

Quantiles_Nuc_counts = quantile(FittedOcc, [0.05:0.05:0.95]);
noQuantiles = size(Quantiles_Nuc_counts, 1);

x = m;
C = [1, .4, .4];
dC = ([1, .85, .85] - C) / 8;

f = figure('Position', [50, 50, 450, 250]);
hold all
for q = [1,5]
    p(q) = patch([x, fliplr(x)],[Quantiles_Nuc_counts(20-q,:), fliplr(Quantiles_Nuc_counts(q,:))], C+(9-q)*dC, 'EdgeColor', 'none');
end
l = plot(x, Quantiles_Nuc_counts(10,:), 'r', 'linewidth', 2);

set(gca, 'FontSize', 11)
xlabel('Digestion time (min)', 'FontSize', 12)
ylabel('Apparent occupancy', 'FontSize', 12)
title('Fitted; 0.3 < O < 0.4', 'FontSize', 12)
h = legend(p([1,5]), {'5%-95%', '25%-75%'}, 'location', 'EO');

v = get(h,'title');
set(v,'string','Percentiles', 'FontSize', 11);

xlim([0 60])
ylim([0 1])
set(gca, 'FontSize', 14)
set(gcf, 'PaperPositionMode', 'auto')
print(gcf, '-depsc', '-painters', fullfile(Figures_folder, 'Percentiles_nuc_counts_low_k1_bin_3_fitted.eps'));


%%
nuc_id = find((k1_vector < 10) & (O_vector >= 0.4) & (O_vector < 0.5));

ratio = norm_factor / max(sum(NucCount, 2));
RescalledNucCount = ratio * NucCount;

Quantiles_Nuc_counts = quantile(RescalledNucCount(nuc_id, :), [0.05:0.05:0.95]);
noQuantiles = size(Quantiles_Nuc_counts, 1);

x = [1, 2, 5, 15, 40, 60];
C = [1, .4, .4];
dC = ([1, .85, .85] - C) / 8;

f = figure('Position', [50, 50, 450, 250]);
hold all
for q = [1,5]
    p(q) = patch([x, fliplr(x)],[Quantiles_Nuc_counts(20-q,:), fliplr(Quantiles_Nuc_counts(q,:))], C+(9-q)*dC, 'EdgeColor', 'none');
end
l = plot(x, Quantiles_Nuc_counts(10,:), 'r', 'linewidth', 2);
    
set(gca, 'FontSize', 11)
xlabel('Digestion time (min)', 'FontSize', 12)
ylabel('Apparent occupancy', 'FontSize', 12)
title('Data; 0.4 < O < 0.5', 'FontSize', 12)
h = legend(p([1,5]), {'5%-95%', '25%-75%'}, 'location', 'EO');

v = get(h,'title');
set(v,'string','Percentiles', 'FontSize', 11);

xlim([0 60])
ylim([0 1])
set(gca, 'FontSize', 14)
set(gcf, 'PaperPositionMode', 'auto')
print(gcf, '-depsc', '-painters', fullfile(Figures_folder, 'Percentiles_nuc_counts_low_k1_bin_4.eps'));

% Plot fitted curves
no_sites = length(nuc_id);
F = @(x, m) x(1) * x(2)/(x(2)-(x(3)))*(exp(-(x(3))*m) - exp(-x(2)*m));

FittedOcc = nan(no_sites, 1001);
m = 0:0.1:100;
for s = 1:no_sites
    FittedOcc(s,:) = F([O_vector(nuc_id(s)), k1_vector(nuc_id(s)), k2_vector(nuc_id(s))], m);
end

Quantiles_Nuc_counts = quantile(FittedOcc, [0.05:0.05:0.95]);
noQuantiles = size(Quantiles_Nuc_counts, 1);

x = m;
C = [1, .4, .4];
dC = ([1, .85, .85] - C) / 8;

f = figure('Position', [50, 50, 450, 250]);
hold all
for q = [1,5]
    p(q) = patch([x, fliplr(x)],[Quantiles_Nuc_counts(20-q,:), fliplr(Quantiles_Nuc_counts(q,:))], C+(9-q)*dC, 'EdgeColor', 'none');
end
l = plot(x, Quantiles_Nuc_counts(10,:), 'r', 'linewidth', 2);

set(gca, 'FontSize', 11)
xlabel('Digestion time (min)', 'FontSize', 12)
ylabel('Apparent occupancy', 'FontSize', 12)
title('Fitted; 0.4 < O < 0.5', 'FontSize', 12)
h = legend(p([1,5]), {'5%-95%', '25%-75%'}, 'location', 'EO');

v = get(h,'title');
set(v,'string','Percentiles', 'FontSize', 11);

xlim([0 60])
ylim([0 1])
set(gca, 'FontSize', 14)
set(gcf, 'PaperPositionMode', 'auto')
print(gcf, '-depsc', '-painters', fullfile(Figures_folder, 'Percentiles_nuc_counts_low_k1_bin_4_fitted.eps'));


%%
nuc_id = find((k1_vector < 10) & (O_vector >= 0.5) & (O_vector < 0.6));

ratio = norm_factor / max(sum(NucCount, 2));
RescalledNucCount = ratio * NucCount;

Quantiles_Nuc_counts = quantile(RescalledNucCount(nuc_id, :), [0.05:0.05:0.95]);
noQuantiles = size(Quantiles_Nuc_counts, 1);

x = [1, 2, 5, 15, 40, 60];
C = [1, .4, .4];
dC = ([1, .85, .85] - C) / 8;

f = figure('Position', [50, 50, 450, 250]);
hold all
for q = [1,5]
    p(q) = patch([x, fliplr(x)],[Quantiles_Nuc_counts(20-q,:), fliplr(Quantiles_Nuc_counts(q,:))], C+(9-q)*dC, 'EdgeColor', 'none');
end
l = plot(x, Quantiles_Nuc_counts(10,:), 'r', 'linewidth', 2);
    
set(gca, 'FontSize', 11)
xlabel('Digestion time (min)', 'FontSize', 12)
ylabel('Apparent occupancy', 'FontSize', 12)
title('Data; 0.5 < O < 0.6', 'FontSize', 12)
h = legend(p([1,5]), {'5%-95%', '25%-75%'}, 'location', 'EO');

v = get(h,'title');
set(v,'string','Percentiles', 'FontSize', 11);

xlim([0 60])
ylim([0 1])
set(gca, 'FontSize', 14)
set(gcf, 'PaperPositionMode', 'auto')
print(gcf, '-depsc', '-painters', fullfile(Figures_folder, 'Percentiles_nuc_counts_low_k1_bin_5.eps'));

% Plot fitted curves
no_sites = length(nuc_id);
F = @(x, m) x(1) * x(2)/(x(2)-(x(3)))*(exp(-(x(3))*m) - exp(-x(2)*m));

FittedOcc = nan(no_sites, 1001);
m = 0:0.1:100;
for s = 1:no_sites
    FittedOcc(s,:) = F([O_vector(nuc_id(s)), k1_vector(nuc_id(s)), k2_vector(nuc_id(s))], m);
end

Quantiles_Nuc_counts = quantile(FittedOcc, [0.05:0.05:0.95]);
noQuantiles = size(Quantiles_Nuc_counts, 1);

x = m;
C = [1, .4, .4];
dC = ([1, .85, .85] - C) / 8;

f = figure('Position', [50, 50, 450, 250]);
hold all
for q = [1,5]
    p(q) = patch([x, fliplr(x)],[Quantiles_Nuc_counts(20-q,:), fliplr(Quantiles_Nuc_counts(q,:))], C+(9-q)*dC, 'EdgeColor', 'none');
end
l = plot(x, Quantiles_Nuc_counts(10,:), 'r', 'linewidth', 2);

set(gca, 'FontSize', 11)
xlabel('Digestion time (min)', 'FontSize', 12)
ylabel('Apparent occupancy', 'FontSize', 12)
title('Fitted; 0.5 < O < 0.6', 'FontSize', 12)
h = legend(p([1,5]), {'5%-95%', '25%-75%'}, 'location', 'EO');

v = get(h,'title');
set(v,'string','Percentiles', 'FontSize', 11);

xlim([0 60])
ylim([0 1])
set(gca, 'FontSize', 14)
set(gcf, 'PaperPositionMode', 'auto')
print(gcf, '-depsc', '-painters', fullfile(Figures_folder, 'Percentiles_nuc_counts_low_k1_bin_5_fitted.eps'));


%%
nuc_id = find((k1_vector < 10) & (O_vector >= 0.6) & (O_vector < 0.7));

ratio = norm_factor / max(sum(NucCount, 2));
RescalledNucCount = ratio * NucCount;

Quantiles_Nuc_counts = quantile(RescalledNucCount(nuc_id, :), [0.05:0.05:0.95]);
noQuantiles = size(Quantiles_Nuc_counts, 1);

x = [1, 2, 5, 15, 40, 60];
C = [1, .4, .4];
dC = ([1, .85, .85] - C) / 8;

f = figure('Position', [50, 50, 450, 250]);
hold all
for q = [1,5]
    p(q) = patch([x, fliplr(x)],[Quantiles_Nuc_counts(20-q,:), fliplr(Quantiles_Nuc_counts(q,:))], C+(9-q)*dC, 'EdgeColor', 'none');
end
l = plot(x, Quantiles_Nuc_counts(10,:), 'r', 'linewidth', 2);
    
set(gca, 'FontSize', 11)
xlabel('Digestion time (min)', 'FontSize', 12)
ylabel('Apparent occupancy', 'FontSize', 12)
title('Data; 0.6 < O < 0.7', 'FontSize', 12)
h = legend(p([1,5]), {'5%-95%', '25%-75%'}, 'location', 'EO');

v = get(h,'title');
set(v,'string','Percentiles', 'FontSize', 11);

xlim([0 60])
ylim([0 1])
set(gca, 'FontSize', 14)
set(gcf, 'PaperPositionMode', 'auto')
print(gcf, '-depsc', '-painters', fullfile(Figures_folder, 'Percentiles_nuc_counts_low_k1_bin_6.eps'));

% Plot fitted curves
no_sites = length(nuc_id);
F = @(x, m) x(1) * x(2)/(x(2)-(x(3)))*(exp(-(x(3))*m) - exp(-x(2)*m));

FittedOcc = nan(no_sites, 1001);
m = 0:0.1:100;
for s = 1:no_sites
    FittedOcc(s,:) = F([O_vector(nuc_id(s)), k1_vector(nuc_id(s)), k2_vector(nuc_id(s))], m);
end

Quantiles_Nuc_counts = quantile(FittedOcc, [0.05:0.05:0.95]);
noQuantiles = size(Quantiles_Nuc_counts, 1);

x = m;
C = [1, .4, .4];
dC = ([1, .85, .85] - C) / 8;

f = figure('Position', [50, 50, 450, 250]);
hold all
for q = [1,5]
    p(q) = patch([x, fliplr(x)],[Quantiles_Nuc_counts(20-q,:), fliplr(Quantiles_Nuc_counts(q,:))], C+(9-q)*dC, 'EdgeColor', 'none');
end
l = plot(x, Quantiles_Nuc_counts(10,:), 'r', 'linewidth', 2);

set(gca, 'FontSize', 11)
xlabel('Digestion time (min)', 'FontSize', 12)
ylabel('Apparent occupancy', 'FontSize', 12)
title('Fitted; 0.6 < O < 0.7', 'FontSize', 12)
h = legend(p([1,5]), {'5%-95%', '25%-75%'}, 'location', 'EO');

v = get(h,'title');
set(v,'string','Percentiles', 'FontSize', 11);

xlim([0 60])
ylim([0 1])
set(gca, 'FontSize', 14)
set(gcf, 'PaperPositionMode', 'auto')
print(gcf, '-depsc', '-painters', fullfile(Figures_folder, 'Percentiles_nuc_counts_low_k1_bin_6_fitted.eps'));


%%
nuc_id = find((k1_vector < 10) & (O_vector >= 0.7) & (O_vector < 0.8));

ratio = norm_factor / max(sum(NucCount, 2));
RescalledNucCount = ratio * NucCount;

Quantiles_Nuc_counts = quantile(RescalledNucCount(nuc_id, :), [0.05:0.05:0.95]);
noQuantiles = size(Quantiles_Nuc_counts, 1);

x = [1, 2, 5, 15, 40, 60];
C = [1, .4, .4];
dC = ([1, .85, .85] - C) / 8;

f = figure('Position', [50, 50, 450, 250]);
hold all
for q = [1,5]
    p(q) = patch([x, fliplr(x)],[Quantiles_Nuc_counts(20-q,:), fliplr(Quantiles_Nuc_counts(q,:))], C+(9-q)*dC, 'EdgeColor', 'none');
end
l = plot(x, Quantiles_Nuc_counts(10,:), 'r', 'linewidth', 2);
    
set(gca, 'FontSize', 11)
xlabel('Digestion time (min)', 'FontSize', 12)
ylabel('Apparent occupancy', 'FontSize', 12)
title('Data; 0.7 < O < 0.8', 'FontSize', 12)
h = legend(p([1,5]), {'5%-95%', '25%-75%'}, 'location', 'EO');

v = get(h,'title');
set(v,'string','Percentiles', 'FontSize', 11);

xlim([0 60])
ylim([0 1])
set(gca, 'FontSize', 14)
set(gcf, 'PaperPositionMode', 'auto')
print(gcf, '-depsc', '-painters', fullfile(Figures_folder, 'Percentiles_nuc_counts_low_k1_bin_7.eps'));

% Plot fitted curves
no_sites = length(nuc_id);
F = @(x, m) x(1) * x(2)/(x(2)-(x(3)))*(exp(-(x(3))*m) - exp(-x(2)*m));

FittedOcc = nan(no_sites, 1001);
m = 0:0.1:100;
for s = 1:no_sites
    FittedOcc(s,:) = F([O_vector(nuc_id(s)), k1_vector(nuc_id(s)), k2_vector(nuc_id(s))], m);
end

Quantiles_Nuc_counts = quantile(FittedOcc, [0.05:0.05:0.95]);
noQuantiles = size(Quantiles_Nuc_counts, 1);

x = m;
C = [1, .4, .4];
dC = ([1, .85, .85] - C) / 8;

f = figure('Position', [50, 50, 450, 250]);
hold all
for q = [1,5]
    p(q) = patch([x, fliplr(x)],[Quantiles_Nuc_counts(20-q,:), fliplr(Quantiles_Nuc_counts(q,:))], C+(9-q)*dC, 'EdgeColor', 'none');
end
l = plot(x, Quantiles_Nuc_counts(10,:), 'r', 'linewidth', 2);

set(gca, 'FontSize', 11)
xlabel('Digestion time (min)', 'FontSize', 12)
ylabel('Apparent occupancy', 'FontSize', 12)
title('Fitted; 0.7 < O < 0.8', 'FontSize', 12)
h = legend(p([1,5]), {'5%-95%', '25%-75%'}, 'location', 'EO');

v = get(h,'title');
set(v,'string','Percentiles', 'FontSize', 11);

xlim([0 60])
ylim([0 1])
set(gca, 'FontSize', 14)
set(gcf, 'PaperPositionMode', 'auto')
print(gcf, '-depsc', '-painters', fullfile(Figures_folder, 'Percentiles_nuc_counts_low_k1_bin_7_fitted.eps'));


%%
nuc_id = find((k1_vector < 10) & (O_vector >= 0.8) & (O_vector < 0.9));

ratio = norm_factor / max(sum(NucCount, 2));
RescalledNucCount = ratio * NucCount;

Quantiles_Nuc_counts = quantile(RescalledNucCount(nuc_id, :), [0.05:0.05:0.95]);
noQuantiles = size(Quantiles_Nuc_counts, 1);

x = [1, 2, 5, 15, 40, 60];
C = [1, .4, .4];
dC = ([1, .85, .85] - C) / 8;

f = figure('Position', [50, 50, 450, 250]);
hold all
for q = [1,5]
    p(q) = patch([x, fliplr(x)],[Quantiles_Nuc_counts(20-q,:), fliplr(Quantiles_Nuc_counts(q,:))], C+(9-q)*dC, 'EdgeColor', 'none');
end
l = plot(x, Quantiles_Nuc_counts(10,:), 'r', 'linewidth', 2);
    
set(gca, 'FontSize', 11)
xlabel('Digestion time (min)', 'FontSize', 12)
ylabel('Apparent occupancy', 'FontSize', 12)
title('Data; 0.8 < O < 0.9', 'FontSize', 12)
h = legend(p([1,5]), {'5%-95%', '25%-75%'}, 'location', 'EO');

v = get(h,'title');
set(v,'string','Percentiles', 'FontSize', 11);

xlim([0 60])
ylim([0 1])
set(gca, 'FontSize', 14)
set(gcf, 'PaperPositionMode', 'auto')
print(gcf, '-depsc', '-painters', fullfile(Figures_folder, 'Percentiles_nuc_counts_low_k1_bin_8.eps'));

% Plot fitted curves
no_sites = length(nuc_id);
F = @(x, m) x(1) * x(2)/(x(2)-(x(3)))*(exp(-(x(3))*m) - exp(-x(2)*m));

FittedOcc = nan(no_sites, 1001);
m = 0:0.1:100;
for s = 1:no_sites
    FittedOcc(s,:) = F([O_vector(nuc_id(s)), k1_vector(nuc_id(s)), k2_vector(nuc_id(s))], m);
end

Quantiles_Nuc_counts = quantile(FittedOcc, [0.05:0.05:0.95]);
noQuantiles = size(Quantiles_Nuc_counts, 1);

x = m;
C = [1, .4, .4];
dC = ([1, .85, .85] - C) / 8;

f = figure('Position', [50, 50, 450, 250]);
hold all
for q = [1,5]
    p(q) = patch([x, fliplr(x)],[Quantiles_Nuc_counts(20-q,:), fliplr(Quantiles_Nuc_counts(q,:))], C+(9-q)*dC, 'EdgeColor', 'none');
end
l = plot(x, Quantiles_Nuc_counts(10,:), 'r', 'linewidth', 2);

set(gca, 'FontSize', 11)
xlabel('Digestion time (min)', 'FontSize', 12)
ylabel('Apparent occupancy', 'FontSize', 12)
title('Fitted; 0.8 < O < 0.9', 'FontSize', 12)
h = legend(p([1,5]), {'5%-95%', '25%-75%'}, 'location', 'EO');

v = get(h,'title');
set(v,'string','Percentiles', 'FontSize', 11);

xlim([0 60])
ylim([0 1])
set(gca, 'FontSize', 14)
set(gcf, 'PaperPositionMode', 'auto')
print(gcf, '-depsc', '-painters', fullfile(Figures_folder, 'Percentiles_nuc_counts_low_k1_bin_8_fitted.eps'));


%%
nuc_id = find((k1_vector < 10) & (O_vector >= 0.9));

ratio = norm_factor / max(sum(NucCount, 2));
RescalledNucCount = ratio * NucCount;

Quantiles_Nuc_counts = quantile(RescalledNucCount(nuc_id, :), [0.05:0.05:0.95]);
noQuantiles = size(Quantiles_Nuc_counts, 1);

x = [1, 2, 5, 15, 40, 60];
C = [1, .4, .4];
dC = ([1, .85, .85] - C) / 8;

f = figure('Position', [50, 50, 450, 250]);
hold all
for q = [1,5]
    p(q) = patch([x, fliplr(x)],[Quantiles_Nuc_counts(20-q,:), fliplr(Quantiles_Nuc_counts(q,:))], C+(9-q)*dC, 'EdgeColor', 'none');
end
l = plot(x, Quantiles_Nuc_counts(10,:), 'r', 'linewidth', 2);
    
set(gca, 'FontSize', 11)
xlabel('Digestion time (min)', 'FontSize', 12)
ylabel('Apparent occupancy', 'FontSize', 12)
title('Data; 0.9 < O', 'FontSize', 12)
h = legend(p([1,5]), {'5%-95%', '25%-75%'}, 'location', 'EO');

v = get(h,'title');
set(v,'string','Percentiles', 'FontSize', 11);

xlim([0 60])
ylim([0 1])
set(gca, 'FontSize', 14)
set(gcf, 'PaperPositionMode', 'auto')
print(gcf, '-depsc', '-painters', fullfile(Figures_folder, 'Percentiles_nuc_counts_low_k1_bin_9.eps'));

% Plot fitted curves
no_sites = length(nuc_id);
F = @(x, m) x(1) * x(2)/(x(2)-(x(3)))*(exp(-(x(3))*m) - exp(-x(2)*m));

FittedOcc = nan(no_sites, 1001);
m = 0:0.1:100;
for s = 1:no_sites
    FittedOcc(s,:) = F([O_vector(nuc_id(s)), k1_vector(nuc_id(s)), k2_vector(nuc_id(s))], m);
end

Quantiles_Nuc_counts = quantile(FittedOcc, [0.05:0.05:0.95]);
noQuantiles = size(Quantiles_Nuc_counts, 1);

x = m;
C = [1, .4, .4];
dC = ([1, .85, .85] - C) / 8;

f = figure('Position', [50, 50, 450, 250]);
hold all
for q = [1,5]
    p(q) = patch([x, fliplr(x)],[Quantiles_Nuc_counts(20-q,:), fliplr(Quantiles_Nuc_counts(q,:))], C+(9-q)*dC, 'EdgeColor', 'none');
end
l = plot(x, Quantiles_Nuc_counts(10,:), 'r', 'linewidth', 2);

set(gca, 'FontSize', 11)
xlabel('Digestion time (min)', 'FontSize', 12)
ylabel('Apparent occupancy', 'FontSize', 12)
title('Fitted; 0.9 < O', 'FontSize', 12)
h = legend(p([1,5]), {'5%-95%', '25%-75%'}, 'location', 'EO');

v = get(h,'title');
set(v,'string','Percentiles', 'FontSize', 11);

xlim([0 60])
ylim([0 1])
set(gca, 'FontSize', 14)
set(gcf, 'PaperPositionMode', 'auto')
print(gcf, '-depsc', '-painters', fullfile(Figures_folder, 'Percentiles_nuc_counts_low_k1_bin_9_fitted.eps'));
