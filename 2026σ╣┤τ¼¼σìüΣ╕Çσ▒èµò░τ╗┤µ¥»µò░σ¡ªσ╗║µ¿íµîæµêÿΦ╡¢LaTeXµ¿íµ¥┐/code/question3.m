% (*@\textcolor{red}{\textbf{本脚本由参赛队员独立编写。编辑过程中曾使用 DeepSeek}}@*)
% (*@\textcolor{red}{\textbf{（deepseek-v4-flash-thinking，杭州深度求索人工智能基础技术研究有限公司，2026）}}@*)
% (*@\textcolor{red}{\textbf{进行代码风格统一与注释清理，核心算法与公式实现均由参赛队员独立完成并验证。}}@*)
clear; clc;

M = readmatrix('附件3.xlsx', 'Range', 'A3');
t  = M(:,1);
z  = M(:,2);
ac = M(:,3);
I  = M(:,4:19);

idx = ~isnan(t) & ~isnan(z) & ~isnan(ac) & all(~isnan(I),2) & z > 0;
t = t(idx);  z = z(idx);  ac = ac(idx);  I = I(idx,:);

mc = 33000;  mf = 11000;  g = 9.8;  K0 = 0.080;
W = 1001;

z_s = smoothdata(z, 'movmean', W);
zd  = gradient(z_s, t);
zdd = gradient(zd, t);

Fe = (mc + mf)*g + mc*ac - mf*zdd;
Fm = K0 * sum(I .* abs(I), 2) ./ z.^2;

valid = isfinite(Fe) & isfinite(Fm) & abs(Fm) > 1e-6;
eta_all = sum(Fm(valid) .* Fe(valid)) / sum(Fm(valid).^2);

tau = 0.10;
main = valid & t >= t(1) + tau & t <= t(end) - tau;
eta_hat = sum(Fm(main) .* Fe(main)) / sum(Fm(main).^2);

res = Fe(main) - eta_hat * Fm(main);
RMSE  = sqrt(mean(res.^2));
MAE   = mean(abs(res));
R2    = 1 - sum(res.^2) / sum((Fe(main) - mean(Fe(main))).^2);
rRMSE = RMSE / mean(abs(Fe(main))) * 100;

eta_pt = Fe(main) ./ Fm(main);

fprintf('eta_all    = %.6f\n', eta_all);
fprintf('eta_hat    = %.6f\n', eta_hat);
fprintf('eta_mean   = %.6f\n', mean(eta_pt));
fprintf('eta_median = %.6f\n', median(eta_pt));
fprintf('eta_std    = %.6f\n', std(eta_pt));
fprintf('R2    = %.4f\n', R2);
fprintf('RMSE  = %.4f\n', RMSE);
fprintf('MAE   = %.4f\n', MAE);
fprintf('rRMSE = %.4f %%\n', rRMSE);

W_list = [501 1001 2001 5001 10001];
eta_W  = zeros(size(W_list));
for q = 1:numel(W_list)
    zs  = smoothdata(z, 'movmean', W_list(q));
    zd_q  = gradient(zs, t);
    zdd_q = gradient(zd_q, t);
    Fe_q = (mc + mf)*g + mc*ac - mf*zdd_q;
    v = isfinite(Fe_q) & isfinite(Fm) & abs(Fm) > 1e-6;
    eta_W(q) = sum(Fm(v) .* Fe_q(v)) / sum(Fm(v).^2);
end

tau_list = [0 0.05 0.10 0.20 0.50];
eta_tau  = zeros(size(tau_list));
for q = 1:numel(tau_list)
    m = valid & t >= t(1) + tau_list(q) & t <= t(end) - tau_list(q);
    eta_tau(q) = sum(Fm(m) .* Fe(m)) / sum(Fm(m).^2);
end

zU_list = [0.0600 0.0599 0.0595 0.0590 0.0580];
eta_zU  = zeros(size(zU_list));
for q = 1:numel(zU_list)
    m = valid & z <= zU_list(q);
    eta_zU(q) = sum(Fm(m) .* Fe(m)) / sum(Fm(m).^2);
end

disp('  window    eta');
disp([W_list(:)  eta_W(:) ]);
disp('  tau       eta');
disp([tau_list(:) eta_tau(:)]);
disp('  z_upper   eta');
disp([zU_list(:)  eta_zU(:) ]);

if eta_hat < 0.8
    verdict = 'attenuation';
elseif eta_hat > 1.2
    verdict = 'amplification';
else
    verdict = 'normal';
end
fprintf('verdict = %s,  attenuation = %.4f %%\n', verdict, (1 - eta_hat)*100);
