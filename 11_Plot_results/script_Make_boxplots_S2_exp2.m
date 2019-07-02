% Load data
load('../../data/S2_exp2/Fit_results.100-200.NF35.S2_exp2.mat')

Figures_folder = 'Figures';
mkdir(Figures_folder)

% Check which nucleosomes is within a DHS
load('S2_DHS_Kharchenko_et_al_Nature_2011.mat', 'DHS_Filter')
isDHS = zeros(size(Chr));

noNucs = numel(Chr);
for n = 1:noNucs
    isDHS(n) = DHS_Filter{Chr(n)}(Loc(n));
end

%% k1 boxplots
figure
data2plot = k1_vector;
h = boxplot(data2plot, isDHS, 'Notch', 'on', 'Whisker', 1, 'OutlierSize', 1);
set(h(7,:), 'Visible', 'off');
set(h,{'linew'},{2})

ylim([0 3.5])
title('k_1 distribution')
xlabel('DHS')
set(gca, 'xticklabel', {'No', 'Yes'});
ylabel('k_1')
set(gca, 'FontSize', 20)

set(gcf, 'PaperPositionMode', 'auto')
print(gcf, '-depsc', '-painters', fullfile(Figures_folder, 'k1_boxplot_DHS.S2_exp2.eps'));
 
%% O boxplots
figure
data2plot = O_vector;
h = boxplot(data2plot, isDHS, 'Notch', 'on', 'Whisker', 1, 'outlierSize', 1);
set(h(7,:), 'Visible', 'off');
set(h,{'linew'},{2})

ylim([0 1])
title('Occupancy distribution')
xlabel('DHS')
set(gca, 'xticklabel', {'No', 'Yes'});
ylabel('Occupancy')
set(gca, 'FontSize', 20)

set(gcf, 'PaperPositionMode', 'auto')
print(gcf, '-depsc', '-painters', fullfile(Figures_folder, 'O_boxplot_DHS.S2_exp2.eps'));

%% k2 boxplots
figure
data2plot = k2_vector;
h = boxplot(data2plot, isDHS, 'Notch', 'on', 'Whisker', 1, 'outlierSize', 1);
set(h(7,:), 'Visible', 'off');
set(h,{'linew'},{2})

ylim([0 0.09])
title('k_2 distribution')
xlabel('DHS')
set(gca, 'xticklabel', {'No', 'Yes'});
ylabel('k_2')
set(gca, 'FontSize', 20)

set(gcf, 'PaperPositionMode', 'auto')
print(gcf, '-depsc', '-painters', fullfile(Figures_folder, 'k2_boxplot_DHS.S2_exp2.eps'));


%%%%%%%%%%%%%%%%%%%%%%%%
%% 5 chromatin states %%
%%%%%%%%%%%%%%%%%%%%%%%%
% Load state annotations
load('Five_Chromatin_States_Filter_dm6.mat', 'ChrStateFilter');
noNucs = numel(Chr);
State = nan(size(Chr));
for n = 1:noNucs
    State(n) = ChrStateFilter{Chr(n)}(Loc(n));
end

noNucs = numel(Chr);

%% k1 boxplots
figure
data2plot = k1_vector;

hold on
h = boxplot(data2plot, State, 'Notch', 'on', 'Whisker', 1, 'OutlierSize', 1);
set(h(7,:), 'Visible', 'off');
set(h,{'linew'},{2})

ylim([0 1.65])
title('k_1 distribution')
xlabel('State')
set(gca, 'xticklabel', {'GREEN', 'YELLOW', 'RED', 'BLUE', 'BLACK'})
xtickangle(45)
ylabel('k_1')
set(gca, 'FontSize', 20)

set(gcf, 'PaperPositionMode', 'auto')
print(gcf, '-depsc', '-painters', fullfile(Figures_folder, 'k1_boxplot_5_states.S2_exp2.eps'));

%% O boxplots
figure
data2plot = O_vector;
h = boxplot(data2plot, State, 'Notch', 'on', 'Whisker', 1, 'outlierSize', 1);
set(h(7,:), 'Visible', 'off');
set(h,{'linew'},{2})

ylim([0 1])
title('Occupancy distribution')
xlabel('State')
set(gca, 'xticklabel', {'GREEN', 'YELLOW', 'RED', 'BLUE', 'BLACK'})
xtickangle(45)
ylabel('Occupancy')
set(gca, 'FontSize', 20)

set(gcf, 'PaperPositionMode', 'auto')
print(gcf, '-depsc', '-painters', fullfile(Figures_folder, 'O_boxplot_5_states.S2_exp2.eps'));

%% k2 boxplots
figure
data2plot = k2_vector;
h = boxplot(data2plot, State, 'Notch', 'on', 'Whisker', 1, 'outlierSize', 1);
set(h(7,:), 'Visible', 'off');
set(h,{'linew'},{2})

ylim([0 0.038])
title('k_2 distribution')
xlabel('State')
set(gca, 'xticklabel', {'GREEN', 'YELLOW', 'RED', 'BLUE', 'BLACK'})
xtickangle(45)
ylabel('k_2')
set(gca, 'FontSize', 20)

set(gcf, 'PaperPositionMode', 'auto')
print(gcf, '-depsc', '-painters', fullfile(Figures_folder, 'k2_boxplot_5_states.S2_exp2.eps'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Promoter gene tertiles %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load promoter annotations
load('PromoterFilter_tx_tertiles.mat', 'PromoterFilter')

noNucs = numel(Chr);
State = nan(size(Chr));
for n = 1:noNucs
    State(n) = PromoterFilter{Chr(n)}(Loc(n));
end

%% k1 boxplots
figure
data2plot = k1_vector;

hold on
h = boxplot(data2plot, State, 'Notch', 'on', 'Whisker', 1, 'OutlierSize', 1);
set(h(7,:), 'Visible', 'off');
set(h,{'linew'},{2})

ylim([0 2.75])
title('k_1 distribution')
set(gca, 'xtick', [1, 2, 3], 'xticklabel', {'High', 'Medium', 'Low'});
xlabel('Transcription level');

% xtickangle(45)
ylabel('k_1')
set(gca, 'FontSize', 20)

set(gcf, 'PaperPositionMode', 'auto')
print(gcf, '-depsc', '-painters', fullfile(Figures_folder, 'k1_boxplot_promoter_3_tertiles.S2_exp2.eps'));

%% O boxplots
figure
data2plot = O_vector;

hold on
h = boxplot(data2plot, State, 'Notch', 'on', 'Whisker', 1, 'OutlierSize', 1);
set(h(7,:), 'Visible', 'off');
set(h,{'linew'},{2})

ylim([0 1])
title('Occupancy distribution')
set(gca, 'xtick', [1, 2, 3], 'xticklabel', {'High', 'Medium', 'Low'});
xlabel('Transcription level');

% xtickangle(45)
ylabel('Occupancy')
set(gca, 'FontSize', 20)

set(gcf, 'PaperPositionMode', 'auto')
print(gcf, '-depsc', '-painters', fullfile(Figures_folder, 'O_boxplot_promoter_3_tertiles.S2_exp2.eps'));

%% k2 boxplots
figure
data2plot = k2_vector;

hold on
h = boxplot(data2plot, State, 'Notch', 'on', 'Whisker', 1, 'OutlierSize', 1);
set(h(7,:), 'Visible', 'off');
set(h,{'linew'},{2})

ylim([0 0.12])
title('k_2 distribution')
set(gca, 'xtick', [1, 2, 3], 'xticklabel', {'High', 'Medium', 'Low'});
xlabel('Transcription level');

% xtickangle(45)
ylabel('k_2')
set(gca, 'FontSize', 20)

set(gcf, 'PaperPositionMode', 'auto')
print(gcf, '-depsc', '-painters', fullfile(Figures_folder, 'k2_boxplot_promoter_3_tertiles.S2_exp2.eps'));


