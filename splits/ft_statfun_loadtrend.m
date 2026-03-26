function [s, cfg] = ft_statfun_loadtrend(cfg, dat, design)

% required
cfg.ivar = ft_getopt(cfg, 'ivar', 1);
cfg.uvar = ft_getopt(cfg, 'uvar', 2);

% find condition indices
cond1 = 1; % load 2
cond2 = 2; % load 4
cond3 = 3; % load 6

sel1 = find(design(cfg.ivar,:) == cond1);
sel2 = find(design(cfg.ivar,:) == cond2);
sel3 = find(design(cfg.ivar,:) == cond3);

% subject IDs
subj = unique(design(cfg.uvar,:));
nsubj = length(subj);

% allocate
contrast = nan(size(dat,1), nsubj);

% compute subject-wise contrast
for i = 1:nsubj
    sidx = subj(i);
    
    idx1 = sel1(design(cfg.uvar,sel1)==sidx);
    idx2 = sel2(design(cfg.uvar,sel2)==sidx);
    idx3 = sel3(design(cfg.uvar,sel3)==sidx);
    
    if isempty(idx1) || isempty(idx2) || isempty(idx3)
        continue
    end
    % average within subject & condition
    d1 = mean(dat(:,idx1),2,'omitnan');
    d2 = mean(dat(:,idx2),2,'omitnan');
    d3 = mean(dat(:,idx3),2,'omitnan');
    
    % linear trend contrast [-1 0 1]
    contrast(:,i) = -d1 + d3;
end

% compute t-statistic across subjects
mean_c = mean(contrast,2,'omitnan');
std_c  = std(contrast,0,2,'omitnan');
n_eff  = sum(isfinite(contrast),2);
denom  = std_c ./ sqrt(max(n_eff, 1));
s.stat = mean_c ./ denom;
s.stat(~isfinite(s.stat)) = 0;

% scalar degrees of freedom expected by FieldTrip statfun interface
s.df = max(nsubj - 1, 1);

% critical value (needed for clustering)
if isfield(cfg, 'alpha')
    s.critval = tinv(1 - cfg.alpha/2, s.df); % two-sided
end

end