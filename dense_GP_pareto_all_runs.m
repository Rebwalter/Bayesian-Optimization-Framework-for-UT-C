%% DENSE GP-BASED PARETO FRONT, ALL 10 RUNS
%
% For each of the 10 stability runs (run_parego_BO_10runs.m): trains
% GP surrogates on that run's own archive, densely samples its Pareto
% front via Chebyshev scalarization + multi-start fmincon (10000
% weight vectors), filters to non-dominated, then computes hypervolume
% against the SAME fixed reference point used everywhere in the paper,
% so it is directly comparable to the raw archive-front hypervolume
% and to the LHS baseline's hypervolume.
%
% Requires: +data/+results/Initial_LHS_Dataset.mat,
%           +data/+results/MOBO_Dataset.mat (from
%           run_parego_BO_10runs.m), helper.paretoFront.m, helper.hypervolume3D.m,
%           helper.area2D_dominated.m, helper.chebyshevCost.m, getFixedRefPoint.m,
%           Parallel Computing Toolbox, Optimization Toolbox,
%           Statistics and Machine Learning Toolbox
%
% Outputs:
%   +data/+results/GP_based_ParetoFronts_Dataset.mat
%   figures/pareto_fronts_10runs_sigma.pdf / .png

clear; clc;

ref_point_fixed = getFixedRefPoint();   % same reference point used for every run/front

load('+data/+results/Initial_LHS_Dataset.mat')
load('+data/+results/MOBO_Dataset.mat')

n_runs = numel(results);

if isempty(gcp('nocreate'))
    parpool('local');
end
poolObj  = gcp();
nWorkers = poolObj.NumWorkers - 1;   % leave one worker free for the client-side loop over runs

paramNames   = fieldnames(Params);
gpParamNames = paramNames(~strcmp(paramNames, 'Radius_tree'));
nGP = numel(gpParamNames);

lb_global = zeros(1,nGP);
ub_global = zeros(1,nGP);
for j = 1:nGP
    p = gpParamNames{j};
    lb_global(j) = Params.(p).lb;
    ub_global(j) = Params.(p).ub;
end

n_weights   = 10000;
n_starts    = 20;
n_starts_ws = 8;
minDist     = 0.10;
rho         = 0.1;   % densification-stage value; see note in runParEGO_BO.m

opts = optimoptions('fmincon', ...
    'Display','off','Algorithm','sqp', ...
    'MaxFunctionEvaluations',500, ...
    'OptimalityTolerance',1e-8, ...
    'StepTolerance',1e-10, ...
    'ConstraintTolerance',1e-8);

%% ----- storage across runs -----
GP_dense_results = struct( ...
    'Y_cheby_pareto',      cell(n_runs,1), ...
    'X_cheby_pareto',      cell(n_runs,1), ...
    'sigmaTotal',          cell(n_runs,1), ...
    'HV_dense_fixedref',   cell(n_runs,1), ...
    'HV_archive_fixedref', cell(n_runs,1), ...
    'n_pareto',            cell(n_runs,1));

%% ---- Sobol weights, sorted for warm-start locality, chunked for parfor ----
% Reused as-is across all 10 runs (the weight vectors themselves don't
% depend on the run; only the GP surrogates being scalarized do).
p = sobolset(2);
p = scramble(p,'MatousekAffineOwen');
sobolWeights = net(p,n_weights);
[sobolWeights_sorted, sortWeightIdx] = sortrows(sobolWeights, [1 2]);

for r = 1:n_runs

    fprintf('\n===== Run %d / %d: training GPs + dense Chebyshev sampling =====\n', r, n_runs);

    %% ---- Build this run's design table + objectives ----
    DOE_r   = results{r}.DOE;
    y_all_r = results{r}.y_all;      % [y_night, y_peak, UTCI_peak]

    X_r = struct2table(DOE_r);
    if any(strcmp(X_r.Properties.VariableNames,'Radius_tree'))
        X_r.Radius_tree = [];
    end

    y_night_r = y_all_r(:,1);
    y_peak_r  = y_all_r(:,2);
    y_utci_r  = y_all_r(:,3);

    %% ---- Train GP surrogates on this run's archive ----
    Mdl_day   = fitrgp(X_r, y_peak_r,  'KernelFunction','ardsquaredexponential','Standardize',true);
    Mdl_night = fitrgp(X_r, y_night_r, 'KernelFunction','ardsquaredexponential','Standardize',true);
    Mdl_utci  = fitrgp(X_r, y_utci_r,  'KernelFunction','ardsquaredexponential','Standardize',true);

    %% ---- Candidate cloud + normalization (per run, since y_all differs) ----
    X_cand_unit = lhsdesign(10000,nGP);
    X_cand = lb_global + (ub_global-lb_global).*X_cand_unit;

    y_min  = min(y_all_r,[],1);
    y_max  = max(y_all_r,[],1);
    y_span = max(y_max-y_min,1e-9);

    yNight = predict(Mdl_night,X_cand);
    yDay   = predict(Mdl_day,X_cand);
    yUTCI  = predict(Mdl_utci,X_cand);
    Ynorm  = ([yNight yDay yUTCI]-y_min)./y_span;

    chunkEdges = round(linspace(0, n_weights, nWorkers+1));

    X_cheby_chunks = cell(nWorkers,1);
    Y_cheby_chunks = cell(nWorkers,1);

    parfor c = 1:nWorkers

        idxRange = (chunkEdges(c)+1):chunkEdges(c+1);
        nLocal = numel(idxRange);

        X_local = zeros(nLocal, nGP);
        Y_local = zeros(nLocal, 3);

        prev_x = [];

        for kk = 1:nLocal

            uvec = sobolWeights_sorted(idxRange(kk), :);
            u1 = uvec(1); u2 = uvec(2);
            lambda = [1-sqrt(u1), sqrt(u1)*(1-u2), sqrt(u1)*u2];

            J = max(Ynorm.*lambda,[],2) + rho*sum(Ynorm.*lambda,2);
            [~,sortIdx] = sort(J,'ascend');

            if isempty(prev_x)
                nFreshNeeded = n_starts;
            else
                nFreshNeeded = n_starts_ws;
            end

            Xstart = [];
            Xstart_norm = [];
            for k = 1:length(sortIdx)
                x = X_cand(sortIdx(k),:);
                x_norm = (x-lb_global)./(ub_global-lb_global);
                if isempty(Xstart)
                    Xstart = x; Xstart_norm = x_norm;
                else
                    d = vecnorm(Xstart_norm-x_norm,2,2);
                    if all(d > minDist)
                        Xstart = [Xstart; x];
                        Xstart_norm = [Xstart_norm; x_norm];
                    end
                end
                if size(Xstart,1) == nFreshNeeded
                    break
                end
            end

            if ~isempty(prev_x)
                Xstart = [prev_x; Xstart];
            end

            objfun = @(x) helper.chebyshevCost(x, Mdl_night, Mdl_day, Mdl_utci, y_min, y_span, lambda, rho);

            best_f = Inf; best_x = []; no_improve = 0;
            for k = 1:size(Xstart,1)
                [x_opt,f_opt] = fmincon(objfun, Xstart(k,:), [],[],[],[], ...
                    lb_global, ub_global, [], opts);
                if f_opt < best_f - 1e-8
                    best_f = f_opt; best_x = x_opt; no_improve = 0;
                else
                    no_improve = no_improve + 1;
                end
                if no_improve >= 5
                    break
                end
            end

            X_local(kk,:) = best_x;
            Y_local(kk,:) = [predict(Mdl_night,best_x), predict(Mdl_day,best_x), predict(Mdl_utci,best_x)];
            prev_x = best_x;
        end

        X_cheby_chunks{c} = X_local;
        Y_cheby_chunks{c} = Y_local;
    end

    X_cheby_sorted = cell2mat(X_cheby_chunks);
    Y_cheby_sorted = cell2mat(Y_cheby_chunks);

    X_cheby = zeros(n_weights, nGP);
    Y_cheby = zeros(n_weights, 3);
    X_cheby(sortWeightIdx,:) = X_cheby_sorted;
    Y_cheby(sortWeightIdx,:) = Y_cheby_sorted;

    %% ---- Filter to non-dominated ----
    mask = helper.paretoFront(Y_cheby);
    Y_cheby_pareto = Y_cheby(mask,:);
    X_cheby_pareto = X_cheby(mask,:);

    %% ---- GP predictive uncertainty on this run's dense Pareto front ----
    [~,sigmaNight_r] = predict(Mdl_night, X_cheby_pareto);
    [~,sigmaDay_r]   = predict(Mdl_day,   X_cheby_pareto);
    [~,sigmaUTCI_r]  = predict(Mdl_utci,  X_cheby_pareto);
    sigmaTotal_r = sqrt(sigmaNight_r.^2 + sigmaDay_r.^2 + sigmaUTCI_r.^2);

    %% ---- Hypervolume: dense GP front vs. raw archive front, SAME fixed ref point ----
    HV_dense = helper.hypervolume3D(Y_cheby_pareto, ref_point_fixed);

    archive_front = y_all_r(results{r}.pareto_idx,:);
    HV_archive = helper.hypervolume3D(archive_front, ref_point_fixed);

    fprintf('Run %d | dense GP front: %d non-dominated / %d | HV=%.4f | archive HV=%.4f\n', ...
        r, sum(mask), n_weights, HV_dense, HV_archive);

    GP_dense_results(r).Y_cheby_pareto      = Y_cheby_pareto;
    GP_dense_results(r).X_cheby_pareto      = X_cheby_pareto;
    GP_dense_results(r).sigmaTotal          = sigmaTotal_r;
    GP_dense_results(r).HV_dense_fixedref   = HV_dense;
    GP_dense_results(r).HV_archive_fixedref = HV_archive;
    GP_dense_results(r).n_pareto            = sum(mask);

end

%% ===================== CROSS-RUN SUMMARY =====================
HV_dense_all   = [GP_dense_results.HV_dense_fixedref]';
HV_archive_all = [GP_dense_results.HV_archive_fixedref]';

fprintf('\n--- Dense GP-front HV (fixed ref) across %d runs ---\n', n_runs);
fprintf('mean=%.4f, std=%.4f, CV=%.2f%%\n', ...
    mean(HV_dense_all), std(HV_dense_all), 100*std(HV_dense_all)/mean(HV_dense_all));

fprintf('\n--- Raw archive-front HV (fixed ref) across %d runs ---\n', n_runs);
fprintf('mean=%.4f, std=%.4f, CV=%.2f%%\n', ...
    mean(HV_archive_all), std(HV_archive_all), 100*std(HV_archive_all)/mean(HV_archive_all));

if ~exist('+data/+results', 'dir')
    mkdir('+data/+results');
end
% save('+data/+results/GP_based_ParetoFronts_Dataset.mat', ...
%     'GP_dense_results', 'ref_point_fixed', '-v7.3')

%% ===================== 3D PARETO FRONTS ACROSS ALL RUNS, COLORED BY SIGMA =====================
% Per-panel titles are kept here (unlike single-panel paper figures,
% which carry no in-figure title) because they are the only thing
% identifying which run each panel belongs to.
fs = 10;

% ---- shared color scale across all panels: global min/max sigma ----
allSigma = vertcat(GP_dense_results.sigmaTotal);
cLimShared = [min(allSigma), max(allSigma)];

% ---- shared axis limits across all panels, for visual comparability ----
% Set from the observed range of this dataset; re-check against your
% own y_all range before reusing with different data.
xLimShared = [13.5, 17.5];
yLimShared = [23.5, 28];
zLimShared = [26, 34];

nCols = 4;
nRows = ceil(n_runs/nCols);

fig = figure('Units','centimeters','Position',[0 0 30 7*nRows],'Color','w');
tiledlayout(nRows,nCols,'TileSpacing','compact','Padding','compact');

for r = 1:n_runs
    nexttile
    scatter3(GP_dense_results(r).Y_cheby_pareto(:,1), ...
             GP_dense_results(r).Y_cheby_pareto(:,2), ...
             GP_dense_results(r).Y_cheby_pareto(:,3), ...
             25, GP_dense_results(r).sigmaTotal, 'filled')

    clim(cLimShared)
    colormap(parula)

    xlabel('T_{night}',  'FontName','Times New Roman','FontSize',fs,'Interpreter','tex')
    ylabel('T_{peak}',   'FontName','Times New Roman','FontSize',fs,'Interpreter','tex')
    zlabel('UTCI_{peak}','FontName','Times New Roman','FontSize',fs,'Interpreter','tex')

    xlim(xLimShared);
    ylim(yLimShared);
    zlim(zLimShared);

    title(sprintf('Run %d (n=%d)', r, GP_dense_results(r).n_pareto), ...
        'FontName','Times New Roman','FontSize',fs,'FontWeight','normal')

    ax = gca;
    ax.FontName  = 'Times New Roman';
    ax.FontSize  = fs-1;
    ax.Box       = 'on';
    ax.LineWidth = 0.75;
    ax.TickDir   = 'out';
    grid on
    ax.GridAlpha = 0.15;
    view(135,25)
end

cb = colorbar;
cb.Layout.Tile    = 'east';
cb.Label.String   = 'GP predictive \sigma (total)';
cb.Label.FontName = 'Times New Roman';
cb.Label.FontSize = fs+1;
cb.FontName       = 'Times New Roman';

if ~exist('figures', 'dir')
    mkdir('figures');
end
exportgraphics(fig, '+figures/pareto_fronts_10runs_sigma.png', 'Resolution',600);
