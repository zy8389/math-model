clear; clc;

M = readmatrix('附件4.xlsx', 'Range', 'A3');
M = M(all(~isnan(M(:,1:19)), 2), 1:19);

t  = M(:,1);
z  = M(:,2);
ac = M(:,3);
I  = M(:,4:19);

idx = z > 0;
t = t(idx);  z = z(idx);  ac = ac(idx);  I = I(idx,:);
[t, ia] = unique(t, 'stable');
z  = z(ia);  ac = ac(ia);  I = I(ia,:);

mc = 33000;  mf = 11000;  g = 9.8;  K0 = 0.080;
W = 1001;
N_mag = 16;
minDur   = 0.50;
mergeGap = 0.10;
tIgn     = 0.50;

thrD =  log(1/0.8);
thrS = -log(1.2);

z_s = smoothdata(z, 'movmean', W);
zd  = gradient(z_s, t);
zdd = gradient(zd, t);

Fe  = (mc + mf)*g + mc*ac - mf*zdd;
Fmi = K0 .* I .* abs(I) ./ z.^2;
Fm  = sum(Fmi, 2);

rho   = Fe ./ Fm;
rho_s = smoothdata(rho, 'movmean', W);
vmid  = t >= t(1) + tIgn & t <= t(end) - tIgn & isfinite(rho) & isfinite(rho_s);

fprintf('rho  mean=%.4f  median=%.4f  min=%.4f  max=%.4f  out=%.4f\n', ...
    mean(rho(vmid)),   median(rho(vmid)), ...
    min(rho(vmid)),    max(rho(vmid)),    ...
    mean(rho(vmid) < 0.8 | rho(vmid) > 1.2));

r  = I.^2 ./ z.^2;
lr = log(r);

e = zeros(size(lr));
for i = 1:N_mag
    e(:,i) = lr(:,i) - median(lr(:, setdiff(1:N_mag, i)), 2);
end
e_s    = smoothdata(e, 1, 'movmean', W);
eta_eq = exp(-e_s);

vdet = t >= t(1) + tIgn & t <= t(end) - tIgn;

rows = {};
for i = 1:N_mag
    flagD = (e_s(:,i) >  thrD) & vdet;
    segD  = findSeg(flagD, t, minDur, mergeGap);
    for k = 1:size(segD, 1)
        ii = segD(k,3):segD(k,4);
        rows(end+1,:) = {i, segD(k,1), segD(k,2), segD(k,2)-segD(k,1), ...
            mean(eta_eq(ii,i)), median(eta_eq(ii,i)), 'decay'};
    end

    flagS = (e_s(:,i) <  thrS) & vdet;
    segS  = findSeg(flagS, t, minDur, mergeGap);
    for k = 1:size(segS, 1)
        ii = segS(k,3):segS(k,4);
        rows(end+1,:) = {i, segS(k,1), segS(k,2), segS(k,2)-segS(k,1), ...
            mean(eta_eq(ii,i)), median(eta_eq(ii,i)), 'surge'};
    end
end

if isempty(rows)
    faultTable = table([],[],[],[],[],[],cell(0,1), ...
        'VariableNames', {'id','t_start','t_end','dur','eta_mean','eta_med','type'});
else
    faultTable = cell2table(rows, 'VariableNames', ...
        {'id','t_start','t_end','dur','eta_mean','eta_med','type'});
    faultTable = sortrows(faultTable, {'t_start','id'});
end
disp(faultTable);

scales = [1.0 1.1 1.2];
sens   = {};
for q = 1:numel(scales)
    s = scales(q);
    ids = [];  cnt = 0;
    for i = 1:N_mag
        flag = ((e_s(:,i) > thrD*s) | (e_s(:,i) < thrS*s)) & vdet;
        seg  = findSeg(flag, t, minDur, mergeGap);
        if ~isempty(seg)
            cnt = cnt + size(seg, 1);
            ids = [ids repmat(i, 1, size(seg,1))];
        end
    end
    sens(end+1,:) = {s, thrD*s, thrS*s, cnt, mat2str(unique(ids))};
end
sensTable = cell2table(sens, 'VariableNames', ...
    {'scale','thr_decay','thr_surge','count','ids'});
disp(sensTable);


function seg = findSeg(flag, t, minDur, mergeGap)
    flag = flag(:);  t = t(:);
    d  = diff([false; flag; false]);
    s0 = find(d ==  1);
    e0 = find(d == -1) - 1;

    seg = [];
    for k = 1:numel(s0)
        if t(e0(k)) - t(s0(k)) >= minDur
            seg = [seg; t(s0(k)), t(e0(k)), s0(k), e0(k)];
        end
    end
    if isempty(seg), return; end

    merged = seg(1,:);
    for k = 2:size(seg, 1)
        if seg(k,1) - merged(end,2) <= mergeGap
            merged(end,2) = seg(k,2);
            merged(end,4) = seg(k,4);
        else
            merged = [merged; seg(k,:)];
        end
    end
    seg = merged;
end
