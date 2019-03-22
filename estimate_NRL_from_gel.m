% Using GelAnalyzer (http://www.gelanalyzer.com/) we obtained the following
% bands corresponding to Sample 3: 191 bp, 386 bp, 579 bp, 791 bp, 1014 bp

% We do a linear regression to find the nucleosome repeat length
bands = [1,2,3,4,5];
lengths = [191, 386, 579, 791, 1014];
fitlm(bands, lengths, 'Intercept', true)

% Estimated Coefficients:
%                    Estimate      SE       tStat       pValue  
%                    ________    ______    _______    __________
% 
%     (Intercept)    -23.1       12.567    -1.8382       0.16333
%     x1             205.1        3.789      54.13    1.3887e-05
% 
% 
% Number of observations: 5, Error degrees of freedom: 3
% Root Mean Squared Error: 12
% R-squared: 0.999,  Adjusted R-Squared 0.999
% F-statistic vs. constant model: 2.93e+03, p-value = 1.39e-05

% We obtained an estimate NRL of 205 bp 
% Standard error of the estimate: 4 bp
