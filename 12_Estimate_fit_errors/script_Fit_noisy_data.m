%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check a single example %
%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate a random example
figure('Position', [50, 50, 700, 300])
hold all
set(gca, 'FontSize', 14)

% Create an example
Occ = 0.8;
K1 = 0.2;
K2 = 0.01;
m = [1:100];
y = Occ * K1/(K1-K2)*(exp(-K2*m) - exp(-K1*m));

% Add noise - uniform distribution between -y/10 and +y/10 (10% relative error)
% Set a seed for the random number generation
seed = 101;
rng(seed);
y_noise = y .* (0.9 + 0.2*rand(1, length(y)));

l(1) = plot(m, y, 'k-', 'linewidth', 2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit using 15 data points %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_points = 15;
sample_points = round(linspace(3,70,N_points));
xdata = m(sample_points);
ydata = y_noise(sample_points);

% Fit parameters: O = x(1); k1 = x(2); k2 = (x(3)/100)
F = @(x, m) x(1) * x(2)/(x(2)-(x(3)/100))*(exp(-(x(3)/100)*m) - exp(-x(2)*m));
options = optimoptions('lsqcurvefit', 'Algorithm', 'trust-region-reflective', ...
    'TolX', 1e-20, 'TolFun', 1e-20, 'MaxFunctionEvaluations', 2000, 'MaxIterations', 1000, 'Display', 'off');
lb = [0.01, 0.01, 0.01];
ub = [1,      50,   50];
x0 = [0.65,  0.6,    1];

x = lsqcurvefit(F, x0, xdata, ydata, lb, ub, options);
y_fit = F(x,m);
l(4) = plot(m, y_fit, 'r', 'linewidth', 1);

plot(m(sample_points), y_noise(sample_points), ...
    'ro', 'MarkerFaceColor','r')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit using 10 data points %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_points = 10;
sample_points = round(linspace(2,65,N_points));
xdata = m(sample_points);
ydata = y_noise(sample_points);

% Fit parameters
x = lsqcurvefit(F, x0, xdata, ydata, lb, ub, options);
y_fit = F(x,m);
l(3) = plot(m, y_fit, 'g', 'linewidth', 1);

plot(m(sample_points), y_noise(sample_points), ...
    'go', 'MarkerFaceColor','g')


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit using 5 data points %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_points = 5;
sample_points = round(linspace(1,60,N_points));
xdata = m(sample_points);
ydata = y_noise(sample_points);

% Fit parameters
x = lsqcurvefit(F, x0, xdata, ydata, lb, ub, options);
y_fit = F(x,m);
l(2) = plot(m, y_fit, 'b', 'linewidth', 1);

plot(m(sample_points), y_noise(sample_points), ...
    'bo', 'MarkerFaceColor','b')


hLegend = legend(l, {'Real function', 'Fitted function (5 data points)', 'Fitted function (10 data points)', 'Fitted function (15 data points)'}, 'location', 'EO');
xlim([0 70])
xlabel('[E]t (arbitrary units)')
ylabel('N/C')

set(gcf, 'PaperPositionMode', 'auto')
print(gcf, '-depsc', '-painters', 'Fit_with_noise.eps');

%% Fit parameters 100000 times and get statistics of the fit
% Set a seed for the random number generation
seed = 101;
rng(seed);

noSims = 100000;

Occ_vector = nan(1, noSims);
K1_vector = nan(1, noSims);
K2_vector = nan(1, noSims);
  
Occ_5 = nan(1, noSims);
K1_5 = nan(1, noSims);
K2_5 = nan(1, noSims);
RMSE_5 = nan(1, noSims);

Occ_10 = nan(1, noSims);
K1_10 = nan(1, noSims);
K2_10 = nan(1, noSims);
RMSE_10 = nan(1, noSims);

Occ_15 = nan(1, noSims);
K1_15 = nan(1, noSims);
K2_15 = nan(1, noSims);
RMSE_15 = nan(1, noSims);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\t Completion: ');
showTimeToCompletion; 
startTime = tic;
p = parfor_progress(noSims);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parfor i = 1:noSims
    
    Occ = 0.5 + 0.5*rand();    % a random value between 10-90 percentiles of the real values
    K1 = 0.3 + 0.7*rand();     % a random value between 10-90 percentiles of the real values
    K2 = 0.002 + 0.027*rand(); % a random value between 10-90 percentiles of the real values
    
    Occ_vector(i) = Occ;
    K1_vector(i) = K1;
    K2_vector(i) = K2;

    m = [1:100];
    y = Occ * K1/(K1-K2)*(exp(-K2*m) - exp(-K1*m));
    y_noise = y .* (0.9 + 0.2*rand(1, length(y)));
    
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Fit using 5 data points %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    sample_points = round(linspace(1,60,5));
    xdata = m(sample_points);
    ydata = y_noise(sample_points);
    
    % Fit parameters
    x = lsqcurvefit(F, x0, xdata, ydata, lb, ub, options);
    yfit = F(x,m)
    RMSE_5(i) = sqrt(mean((y - yfit).^2)); 
    Occ_5(i) = x(1);
    K1_5(i)  = x(2);
    K2_5(i)  = x(3)/100;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Fit using 10 data points %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sample_points = round(linspace(1,60,10));
    xdata = m(sample_points);
    ydata = y_noise(sample_points);
    
    % Fit parameters
    x = lsqcurvefit(F, x0, xdata, ydata, lb, ub, options);
    yfit = F(x,m)
    RMSE_10(i) = sqrt(mean((y - yfit).^2)); 
    Occ_10(i) = x(1);
    K1_10(i)  = x(2);
    K2_10(i)  = x(3)/100;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Fit using 15 data points %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sample_points = round(linspace(1,60,15));
    xdata = m(sample_points);
    ydata = y_noise(sample_points);
    
    % Fit parameters
    x = lsqcurvefit(F, x0, xdata, ydata, lb, ub, options);
    yfit = F(x,m)
    RMSE_15(i) = sqrt(mean((y - yfit).^2)); 
    Occ_15(i) = x(1);
    K1_15(i)  = x(2);
    K2_15(i)  = x(3)/100;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    p = parfor_progress;
    showTimeToCompletion(p/100, [], [], startTime);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%%
Occ_error_15 = abs(Occ_15-Occ_vector)./Occ_vector;
Occ_error_10 = abs(Occ_10-Occ_vector)./Occ_vector;
Occ_error_5 = abs(Occ_5-Occ_vector)./Occ_vector;

K1_error_15  = abs(K1_15-K1_vector)./K1_vector;
K1_error_10  = abs(K1_10-K1_vector)./K1_vector;
K1_error_5  = abs(K1_5-K1_vector)./K1_vector;

K2_error_15  = abs(K2_15-K2_vector)./K2_vector;
K2_error_10  = abs(K2_10-K2_vector)./K2_vector;
K2_error_5  = abs(K2_5-K2_vector)./K2_vector;

save('Sim_results.mat', 'K1_10', 'K1_15', 'K1_5', 'K1_error_10', 'K1_error_15', 'K1_error_5', 'K1_vector', ...
    'K2_10', 'K2_15', 'K2_5', 'K2_error_10', 'K2_error_15', 'K2_error_5', 'K2_vector', ...
    'Occ_10', 'Occ_15', 'Occ_5', 'Occ_error_10', 'Occ_error_15', 'Occ_error_5', 'Occ_vector', ...
    'RMSE_5', 'RMSE_10', 'RMSE_15', 'noSims');


%%
Error = 100*[Occ_error_5, Occ_error_10, Occ_error_15, K1_error_5, K1_error_10, K1_error_15, K2_error_5, K2_error_10, K2_error_15];
Num_of_points = [5*ones(1, noSims), 10*ones(1, noSims), 15*ones(1, noSims), 5*ones(1, noSims), 10*ones(1, noSims), 15*ones(1, noSims), 5*ones(1, noSims), 10*ones(1, noSims), 15*ones(1, noSims)];
[Parameter{[1:3*noSims]}] = deal('O');
[Parameter{3*noSims+[1:3*noSims]}] = deal('K1');
[Parameter{6*noSims+[1:3*noSims]}] = deal('K2');

clear g
g = gramm('x',Num_of_points, ...
          'y',Error, ...
          'color',Parameter);
          
g.stat_boxplot();
g.axe_property('XLim',[3 17],'YLim',[0 35]);
g.set_names('color','Parameter','x','Number of points','y','Relative error (%)');
figure('Position', [50, 50, 400, 300])
g.draw();

set(gcf, 'PaperPositionMode', 'auto')
print(gcf, '-depsc', '-painters', 'Fit_summary.eps');


%%
median(K1_error_5)  % 0.0897
median(K2_error_5)  % 0.0961
median(Occ_error_5) % 0.0555

median(K1_error_10)  % 0.0763
median(K2_error_10)  % 0.0640
median(Occ_error_10) % 0.0311

median(K1_error_15)  % 0.0765
median(K2_error_15)  % 0.0512
median(Occ_error_15) % 0.0242


%%
corr(Occ_vector(:), Occ_5(:))  % 0.9371
corr(K1_vector(:), K1_5(:))  % 0.9054
corr(K2_vector(:), K2_5(:))  % 0.9712

corr(Occ_vector(:), Occ_10(:)) % 0.9744
corr(K1_vector(:), K1_10(:)) % 0.9378
corr(K2_vector(:), K2_10(:)) % 0.9860

corr(Occ_vector(:), Occ_15(:)) % 0.9833
corr(K1_vector(:), K1_15(:)) % 0.9416
corr(K2_vector(:), K2_15(:)) % 0.9906
