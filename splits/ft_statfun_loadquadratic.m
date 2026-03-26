function [s, cfg] = ft_statfun_loadquadratic(cfg, dat, design)

cfg.ivar = ft_getopt(cfg, 'ivar', 1);
cfg.uvar = ft_getopt(cfg, 'uvar', 2);

cond1 = 1; % load 2
cond2 = 2; % load 4
cond3 = 3; % load 6

sel1 = find(design(cfg.ivar,:) == cond1);
sel2 = find(design(cfg.ivar,:) == cond2);
sel3 = find(design(cfg.ivar,:) == cond3);

subj = unique(design(cfg.uvar,:));
nsubj = length(subj);

contrast = nan(size(dat,1), nsubj);

for i = 1:nsubj
    sidx = subj(i);
    
    idx1 = sel1(design(cfg.uvar,sel1)==sidx);
    idx2 = sel2(design(cfg.uvar,sel2)==sidx);
    idx3 = sel3(design(cfg.uvar,sel3)==sidx);
    
    if isempty(idx1) || isempty(idx2) || isempty(idx3)
        continue
    end
    d1 = mean(dat(:,idx1),2,'omitnan');
    d2 = mean(dat(:,idx2),2,'omitnan');
    d3 = mean(dat(:,idx3),2,'omitnan');
    
    % quadratic contrast [1 -2 1]
    contrast(:,i) = d1 - 2*d2 + d3;
end

% t-statistic
mean_c = mean(contrast,2,'omitnan');
std_c  = std(contrast,0,2,'omitnan');
n_eff  = sum(isfinite(contrast),2);
denom  = std_c ./ sqrt(max(n_eff, 1));
s.stat = mean_c ./ denom;
s.stat(~isfinite(s.stat)) = 0;
s.df   = max(nsubj - 1, 1);

if isfield(cfg,'alpha')
    s.critval = tinv(1 - cfg.alpha/2, s.df); % two-sided
end

end