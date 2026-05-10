clear; clc;

M = readmatrix('附件2.xlsx', 'Range', 'A3');
t  = M(:,1);
F16 = M(:,2:17);

idx = ~isnan(t) & all(~isnan(F16),2) & t >= 0 & t <= 10;
t = t(idx);  F16 = F16(idx,:);
[t, ia] = unique(t, 'stable');
F16 = F16(ia,:);

Fe_series = sum(F16, 2);
Fe = @(tt) interp1(t, Fe_series, tt, 'linear', 'extrap');

mc = 33000;  mf = 11000;  g = 9.8;
k  = 2e7;    c  = 8e4;    z_max = 0.06;

N = length(t);
X = zeros(N, 4);
X(1,:) = [-mc*g/k, 0, 0, 0];
ground = false(N, 1);
ground(1) = (mf*g - Fe(t(1)) - k*X(1,1) - c*X(1,2)) >= 0;

for n = 1:N-1
    dt = t(n+1) - t(n);
    xk = X(n,:).';

    if ground(n)
        zc = xk(1);  vc = xk(2);
        R = mf*g - Fe(t(n)) - k*zc - c*vc;

        if R >= 0
            f = @(tt, yy) [yy(2); (-k*yy(1) - c*yy(2) - mc*g) / mc];
            y = rk4(f, t(n), [zc; vc], dt);
            R_next = mf*g - Fe(t(n+1)) - k*y(1) - c*y(2);
            X(n+1,:) = [y(1), y(2), 0, 0];
            ground(n+1) = R_next >= 0;
        else
            f = @(tt, xx) dyn(tt, xx, Fe, mc, mf, k, c, g);
            xn = rk4(f, t(n), [xk(1); xk(2); 0; 0], dt);
            [X(n+1,:), ground(n+1)] = proj(xn, t(n+1), Fe, mf, k, c, g, z_max);
        end
    else
        f = @(tt, xx) dyn(tt, xx, Fe, mc, mf, k, c, g);
        xn = rk4(f, t(n), xk, dt);
        [X(n+1,:), ground(n+1)] = proj(xn, t(n+1), Fe, mf, k, c, g, z_max);
    end
end

zc = X(:,1);  zf = X(:,3);
zc_9  = interp1(t, zc, 9);
zf_9  = interp1(t, zf, 9);
gap_9 = z_max - zf_9;

fprintf('zc(9)    = %.6f m\n', zc_9);
fprintf('zf(9)    = %.6f m\n', zf_9);
fprintf('zgap(9)  = %.6f m  = %.4f mm\n', gap_9, gap_9 * 1000);


function dx = dyn(t, x, Fe, mc, mf, k, c, g)
    zc = x(1);  vc = x(2);
    zf = x(3);  vf = x(4);
    Fs = k*(zf - zc) + c*(vf - vc);
    dx = [vc;
          (Fs - mc*g)/mc;
          vf;
          (Fe(t) - Fs - mf*g)/mf];
end

function y = rk4(f, t, x, dt)
    k1 = f(t,        x);
    k2 = f(t + dt/2, x + dt*k1/2);
    k3 = f(t + dt/2, x + dt*k2/2);
    k4 = f(t + dt,   x + dt*k3);
    y  = x + dt*(k1 + 2*k2 + 2*k3 + k4)/6;
end

function [xp, ground] = proj(x, t, Fe, mf, k, c, g, z_max)
    zc = x(1);  vc = x(2);
    zf = x(3);  vf = x(4);
    ground = false;

    if zf < 0
        zf = 0;
        R = mf*g - Fe(t) - k*zc - c*vc;
        if R >= 0
            ground = true;
            vf = 0;
        else
            if vf < 0, vf = 0; end
        end
    end

    if zf > z_max
        zf = z_max;
        if vf > 0, vf = 0; end
    end

    xp = [zc, vc, zf, vf];
end
