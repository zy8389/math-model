clear; clc;

T = readtable('附件1.xlsx', 'VariableNamingRule', 'preserve');
F = T{:,1};
I = T{:,2};
z = T{:,3};

idx = ~isnan(F) & ~isnan(I) & ~isnan(z);
F = F(idx);  I = I(idx);  z = z(idx);

x = I .* abs(I) ./ z.^2;
K0 = sum(x .* F) / sum(x.^2);

e = F - K0 * x;
R2    = 1 - sum(e.^2) / sum((F - mean(F)).^2);
RMSE  = sqrt(mean(e.^2));
MAE   = mean(abs(e));
rRMSE = RMSE / mean(abs(F)) * 100;

fprintf('K0    = %.8f\n', K0);
fprintf('R2    = %.4f\n', R2);
fprintf('RMSE  = %.4f N\n', RMSE);
fprintf('MAE   = %.4f N\n', MAE);
fprintf('rRMSE = %.4f %%\n', rRMSE);
fprintf('mean_e = %.4e,  std_e = %.4e\n', mean(e), std(e));
fprintf('corr(z,e) = %.4f\n', corr(z, e));
fprintf('corr(I,e) = %.4f\n', corr(I, e));
