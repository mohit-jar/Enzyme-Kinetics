% Enzyme Kinetics Analysis

data = [25 50 100 250 500 1000 2500 5000 7500 10000;  % Substrate concentration
        0.0243 0.0292 0.0546 0.1388 0.1726 0.2374 0.3023 0.3395 0.3515 0.3652];  % Reaction rate

S = data(1, :); % Substrate concentration
V = data(2, :); % Reaction rate

% Define Michaelis-Menten Model
mmModel = @(beta, S) (beta(1) .* S) ./ (beta(2) + S);

% Initial guesses for Vmax and Km
beta0 = [max(V), mean(S)];

% Lineweaver-Burk Transformation
LB_S = 1 ./ S;
LB_V = 1 ./ V;
p_LB = polyfit(LB_S, LB_V, 1);
Km_LB = p_LB(1) / p_LB(2);
Vmax_LB = 1 / p_LB(2);
V_fit_LB = mmModel([Vmax_LB, Km_LB], S);

% Eadie-Hofstee Transformation
EH_V = V;
EH_S = V ./ S;
p_EH = polyfit(EH_S, EH_V, 1);
Vmax_EH = p_EH(2);
Km_EH = -p_EH(1);
V_fit_EH = mmModel([Vmax_EH, Km_EH], S);

% Hanes-Woolf Transformation
HW_S = S;
HW_V = S ./ V;
p_HW = polyfit(HW_S, HW_V, 1);
Km_HW = p_HW(2) / p_HW(1);
Vmax_HW = 1 / p_HW(1);
V_fit_HW = mmModel([Vmax_HW, Km_HW], S);

% Ordinary Least Squares (OLS) using nlinfit without robust option
[beta_ols, R_ols, J_ols, covB_ols, mse_ols] = nlinfit(S, V, mmModel, beta0);

% Robust Nonlinear Regression using nlinfit with bisquare robust option
opts = statset('nlinfit');
opts.RobustWgtFun = 'bisquare';
[beta_robust, R, J, covB, mse] = nlinfit(S, V, mmModel, beta0, opts);

% Calculate 95% confidence intervals for OLS and robust fit parameters
ci_ols = nlparci(beta_ols, R_ols, 'covar', covB_ols);
ci_robust = nlparci(beta_robust, R, 'covar', covB);

% Calculate Mean Absolute Error for each model
MAE_LB = mean(abs(V - V_fit_LB));
MAE_EH = mean(abs(V - V_fit_EH));
MAE_HW = mean(abs(V - V_fit_HW));
MAE_ols = mean(abs(V - mmModel(beta_ols, S)));
MAE_robust = mean(abs(V - mmModel(beta_robust, S)));

% Visualize Results
figure;
subplot(2,2,1);
plot(S, V, 'bo');
hold on;
plot(S, V_fit_LB, 'r-');
title('Lineweaver-Burk');
xlabel('S');
ylabel('V');
legend('Data', 'Lineweaver-Burk Fit');
text(0.5, 0.9, ['MAE = ' num2str(MAE_LB)], 'Units', 'normalized');

subplot(2,2,2);
plot(S, V, 'bo');
hold on;
plot(S, V_fit_EH, 'r-');
title('Eadie-Hofstee');
xlabel('S');
ylabel('V');
legend('Data', 'Eadie-Hofstee Fit');
text(0.5, 0.9, ['MAE = ' num2str(MAE_EH)], 'Units', 'normalized');

subplot(2,2,3);
plot(S, V, 'bo');
hold on;
plot(S, V_fit_HW, 'r-');
title('Hanes-Woolf');
xlabel('S');
ylabel('V');
legend('Data', 'Hanes-Woolf Fit');
text(0.5, 0.9, ['MAE = ' num2str(MAE_HW)], 'Units', 'normalized');

subplot(2,2,4);
plot(S, V, 'bo');
hold on;
plot(S, mmModel(beta0, S), 'g--');
plot(S, mmModel(beta_robust, S), 'k-');
title('Robust Nonlinear Regression');
xlabel('S');
ylabel('V');
legend('Data', 'Initial Guess', 'Robust Fit');
text(0.5, 0.9, ['MAE = ' num2str(MAE_robust)], 'Units', 'normalized');

% Combined Plot (S vs. V with model fits)
figure;
plot(S, V, 'bo');
hold on;
plot(S, V_fit_LB, 'r--');
plot(S, V_fit_EH, 'm--');
plot(S, V_fit_HW, 'c--');
plot(S, mmModel(beta_ols, S), 'b-');
plot(S, mmModel(beta_robust, S), 'k-');
legend('Data', 'Lineweaver-Burk', 'Eadie-Hofstee', 'Hanes-Woolf', 'OLS Fit', 'Robust Fit');
title('Comparison of Different Transformations, OLS, and Robust Nonlinear Regression');
xlabel('S');
ylabel('V');

% Display estimated parameters and confidence intervals
fprintf('Lineweaver-Burk: Km = %.4f, Vmax = %.4f, MAE = %.4f\n', Km_LB, Vmax_LB, MAE_LB);
fprintf('Eadie-Hofstee: Km = %.4f, Vmax = %.4f, MAE = %.4f\n', Km_EH, Vmax_EH, MAE_EH);
fprintf('Hanes-Woolf: Km = %.4f, Vmax = %.4f, MAE = %.4f\n', Km_HW, Vmax_HW, MAE_HW);
fprintf('OLS: Km = %.4f (95%% CI: [%.4f, %.4f]), Vmax = %.4f (95%% CI: [%.4f, %.4f]), MAE = %.4f\n', ...
        beta_ols(2), ci_ols(2,1), ci_ols(2,2), beta_ols(1), ci_ols(1,1), ci_ols(1,2), MAE_ols);
fprintf('Robust Nonlinear Regression: Km = %.4f (95%% CI: [%.4f, %.4f]), Vmax = %.4f (95%% CI: [%.4f, %.4f]), MAE = %.4f\n', ...
        beta_robust(2), ci_robust(2,1), ci_robust(2,2), beta_robust(1), ci_robust(1,1), ci_robust(1,2), MAE_robust);
