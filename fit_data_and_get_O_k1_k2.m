load('data/NucCounts.mat', 'Chr', 'Loc', 'NucCount', 'noNucs');

%% Fit all nuc data and get the occ and acc
O_vector = nan(noNucs, 1);
k1_vector = nan(noNucs, 1);
k2_vector = nan(noNucs, 1);

m = [1, 2, 5, 15, 40, 60];

% In a sequencing experiment we can only obtain the number of nucleosomes
% that were sequenced by the sequencer, and not the real number of
% nucleosomes that were present in the sample, before sequencing, because
% it is impossible to obtain the fraction of DNA fragments that are
% actually captured and read by a sequencer. Therefore, the fraction of
% cells that released a nucleosome, N/C (#Nuc. / #Cells), is not directly
% measurable.
% 
% norm_factor is the proportionality factor that transforms nuc. counts
% into apparent nuc. occupancy, i.e. #Nuc. / #Cells; norm_factor controls
% the average nuc. occupancy. Reasonable values for norm_factor (resulting
% in average nucleosome occupancy 0.4 - 0.8) are: 30-50. For each set of
% sequencing experiments, one needs to test different values of
% norm_factor, which will account for the losses in the specific sequencer
% (the fraction of DNA that is not sequenced). For demonstration purposes,
% below we use a norm_factor of 40.
norm_factor = 40;
ratio = norm_factor / max(sum(NucCount, 2));

parfor c = 1:noNucs
    NucRatio = ratio * NucCount(c,:);
    
    % Fit parameters: O = x(1); k1 = x(2); k2 = (x(3)/100)
    F = @(x, m) x(1) * x(2)/(x(2)-(x(3)/100))*(exp(-(x(3)/100)*m) - exp(-x(2)*m));
    options = optimoptions('lsqcurvefit', 'Algorithm', 'trust-region-reflective', ...
        'TolX', 1e-20, 'TolFun', 1e-20, 'MaxFunctionEvaluations', 2000, 'MaxIterations', 1000, 'Display', 'off');
    xdata = m;
    ydata = NucRatio;
    lb = [0.01, 0.01, 0.01];
    ub = [1,      50,   50];
    x0 = [0.65,  0.6,    1];
    
    x = lsqcurvefit(F, x0, xdata, ydata, lb, ub, options);
    O_vector(c) = x(1);
    k1_vector(c) = x(2);
    k2_vector(c) = x(3)/100;
    
    % disp(c)
end

save('data/Fit_results.mat', 'Chr', 'Loc', 'NucCount', 'O_vector', 'k1_vector', 'k2_vector');
