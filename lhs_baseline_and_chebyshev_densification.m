%% LHS BASELINE + CHEBYSHEV DENSIFICATION (COMPARISON METHOD)
%
% Naive-budget comparison method, evaluated against the ParEGO
% Bayesian-optimization pipeline (generate_initial_LHS_dataset.m ->
% run_parego_BO_10runs.m -> dense_GP_pareto_all_runs.m). This script
% does NOT feed into that pipeline; it is a self-contained baseline
% used to show what a purely-LHS, non-adaptive design achieves for a
% comparable evaluation budget:
%
%   1. Draw a 144-point LHS design and evaluate it with UT&C.
%   2. Compute the non-dominated (Pareto) subset and its hypervolume
%      against the fixed reference point (same one used everywhere).
%   3. Train 3 GP surrogates on the 144 LHS points.
%   4. Densify the GP-approximated Pareto front with 2000
%      Chebyshev-scalarized points (Sobol weights + multi-start
%      fmincon), and check how many of the 144 true (UT&C-evaluated)
%      Pareto points remain non-dominated once compared against that
%      denser GP front.
%
% Requires: helper.makeDOE_opt.m, UTC.m, helper.paretoFront.m,
%           helper.hypervolume3D.m, helper.area2D_dominated.m, helper.chebyshevCost.m,
%           getFixedRefPoint.m, Parallel Computing Toolbox, Statistics
%           and Machine Learning Toolbox (fitrgp, sobolset)
%
% Outputs:
%   +data/+results/LHS_Baseline_Dataset.mat
%   +data/+results/LHS_Baseline_Chebyshev.mat
%   figures/lhs_baseline_pareto_densification.pdf / .png

clear; clc;

%% ===================== SITE / DOE CONFIGURATION =====================
load('+data/+input/MeteoData_DWD_hourly_RadPart.mat');
MeteoData = MeteoData_DWD(1009:1105,:);   % simulation window (hourly)
Name_Site = 'HST_BEM';
m         = 144;   % number of LHS design points
rdmseed   = 11;      % fixed seed for reproducibility

[DOE, ~, Params] = helper.makeDOE_opt(m, rdmseed);

%% ===================== PARALLEL UT&C EVALUATION =====================
if isempty(gcp('nocreate'))
    parpool('local', 11);
end
poolObj  = gcp();
nWorkers = poolObj.NumWorkers;
fprintf('Using %d parallel workers\n', nWorkers);

% Split the DOE into contiguous chunks so worker order still matches the
% original DOE order once reassembled below.
chunkEdges = round(linspace(0, m, nWorkers+1));
DOE_chunks = cell(nWorkers,1);
m_chunks   = zeros(nWorkers,1);
for c = 1:nWorkers
    idxRange      = (chunkEdges(c)+1):chunkEdges(c+1);
    DOE_chunks{c} = DOE(idxRange);
    m_chunks(c)   = numel(idxRange);
end

T2m_chunks  = cell(nWorkers,1);
UTCI_chunks = cell(nWorkers,1);
parfor c = 1:nWorkers
    [Results2m_c, UTCI_c] = UTC(DOE_chunks{c}, m_chunks(c), MeteoData, Name_Site);
    Results2m_c.T2m = Results2m_c.T2m - 273.15;   % K -> degC

    T2m_chunks{c}  = Results2m_c.T2m;
    UTCI_chunks{c} = UTCI_c;
end

Results2m.T2m = cat(3, T2m_chunks{:});
UTCI          = cat(3, UTCI_chunks{:});

assert(size(Results2m.T2m,3) == m, ...
    'Reassembled T2m has %d slices, expected m=%d', size(Results2m.T2m,3), m);
assert(size(UTCI,3) == m, ...
    'Reassembled UTCI has %d slices, expected m=%d', size(UTCI,3), m);

%% ===================== OBJECTIVES =====================
% All three objectives are minimized: hot daytime peak, cold nighttime
% low (a low value = a colder night, the extreme of interest), and
% UTCI peak.
T_peak       = max(Results2m.T2m(:,:,:), [], 1);
y_peak       = reshape(T_peak, 1, [])';

T_window     = Results2m.T2m(40:60,:,:);   % nighttime hours window
NightMinTemp = min(T_window, [], 1);
y_night      = reshape(NightMinTemp, 1, [])';

UTCI_peak    = max(UTCI, [], 1);
UTCI_peak    = reshape(UTCI_peak, 1, [])';

y_all = [y_night, y_peak, UTCI_peak];

%% ===================== NON-DOMINATED (PARETO) SUBSET =====================
isPareto       = helper.paretoFront(y_all);
pareto_idx_LHS = find(isPareto);
Y_pareto_LHS   = y_all(pareto_idx_LHS,:);

fprintf('LHS baseline: %d design points, %d non-dominated points\n', ...
    m, numel(pareto_idx_LHS));

%% ===================== HYPERVOLUME (fixed reference point) =====================
ref_point_fixed = getFixedRefPoint();   % same reference point used for every run/front
HV_LHS = helper.hypervolume3D(Y_pareto_LHS, ref_point_fixed);

fprintf('LHS baseline hypervolume (fixed ref): %.4f\n', HV_LHS);

%% ===================== TRAIN 3 GP SURROGATES =====================
% Radius_tree is a derived quantity (not an independent design
% variable) and is excluded from the surrogate inputs.
paramNames   = fieldnames(Params);
gpParamNames = paramNames(~strcmp(paramNames, 'Radius_tree'));

X = struct2table(DOE);
if any(strcmp(X.Properties.VariableNames, 'Radius_tree'))
    X.Radius_tree = [];
end

y_night_shape = y_all(:,1);
y_peak_shape  = y_all(:,2);
y_UTCI_shape  = y_all(:,3);

Mdl_night = fitrgp(X, y_night_shape, 'KernelFunction','ardsquaredexponential', 'Standardize', true);
Mdl_day   = fitrgp(X, y_peak_shape,  'KernelFunction','ardsquaredexponential', 'Standardize', true);
Mdl_utci  = fitrgp(X, y_UTCI_shape,  'KernelFunction','ardsquaredexponential', 'Standardize', true);

fprintf('GPs trained: Mdl_night, Mdl_day, Mdl_utci (on %d points, %d parameters)\n', ...
    size(X,1), size(X,2));

%% ===================== SAVE LHS BASELINE =====================
if ~exist('+data/+results', 'dir')
    mkdir('+data/+results');
end
% save('+data/+results/LHS_Baseline_Dataset.mat', ...
%     'Params', 'gpParamNames', 'DOE', 'X', 'm', 'Name_Site', 'MeteoData', ...
%     'y_all', 'y_peak', 'y_night', 'UTCI_peak', ...
%     'Y_pareto_LHS', 'pareto_idx_LHS', 'HV_LHS', 'ref_point_fixed', ...
%     'Mdl_night', 'Mdl_day', 'Mdl_utci');

%% ===================== DENSIFY WITH CHEBYSHEV SCALARIZATION =====================
% Uses a fresh worker pool sizing pass since this stage is far more
% compute-heavy per task (multi-start fmincon at every Sobol weight)
% than the UT&C evaluation above.
if isempty(gcp('nocreate'))
    parpool('local');
end
poolObj  = gcp();
nWorkers = poolObj.NumWorkers;

nGP = numel(gpParamNames);
lb_global = zeros(1, nGP);
ub_global = zeros(1, nGP);
for j = 1:nGP
    p = gpParamNames{j};
    lb_global(j) = Params.(p).lb;
    ub_global(j) = Params.(p).ub;
end

n_weights   = 10000;    % number of Sobol weight vectors (Pareto-front density)
n_starts    = 20;      % multi-starts for the first weight in each chunk
n_starts_ws = 8;        % warm-started multi-starts for subsequent weights
minDist     = 0.10;    % minimum normalized distance between multi-start points
rho         = 0.1;     % augmented-Tchebycheff augmentation coefficient (densification-stage value; see note in runParEGO_BO.m)

opts = optimoptions('fmincon', ...
    'Display', 'off', 'Algorithm', 'sqp', ...
    'MaxFunctionEvaluations', 500, ...
    'OptimalityTolerance', 1e-8, ...
    'StepTolerance', 1e-10, ...
    'ConstraintTolerance', 1e-8);

X_cand_unit = lhsdesign(10000, nGP);
X_cand      = lb_global + (ub_global - lb_global) .* X_cand_unit;

y_min  = min(y_all, [], 1);
y_max  = max(y_all, [], 1);
y_span = max(y_max - y_min, 1e-9);

yNight = predict(Mdl_night, X_cand);
yDay   = predict(Mdl_day,   X_cand);
yUTCI  = predict(Mdl_utci,  X_cand);
Ynorm  = ([yNight, yDay, yUTCI] - y_min) ./ y_span;

p = sobolset(2);
p = scramble(p, 'MatousekAffineOwen');
sobolWeights = net(p, n_weights);

% Sort so that adjacent weight vectors are spatially close, which lets
% each chunk warm-start fmincon from the previous weight's solution.
[sobolWeights_sorted, sortWeightIdx] = sortrows(sobolWeights, [1 2]);
chunkEdges = round(linspace(0, n_weights, nWorkers+1));

X_cheby_chunks = cell(nWorkers,1);
Y_cheby_chunks = cell(nWorkers,1);

parfor c = 1:nWorkers
    idxRange = (chunkEdges(c)+1):chunkEdges(c+1);
    nLocal   = numel(idxRange);

    X_local = zeros(nLocal, nGP);
    Y_local = zeros(nLocal, 3);
    prev_x  = [];

    for kk = 1:nLocal
        uvec = sobolWeights_sorted(idxRange(kk), :);
        u1 = uvec(1); u2 = uvec(2);
        % Map the 2D Sobol point to a uniform weight on the 3-simplex.
        lambda = [1-sqrt(u1), sqrt(u1)*(1-u2), sqrt(u1)*u2];

        J = max(Ynorm.*lambda, [], 2) + rho*sum(Ynorm.*lambda, 2);
        [~, sortIdx] = sort(J, 'ascend');

        if isempty(prev_x)
            nFreshNeeded = n_starts;
        else
            nFreshNeeded = n_starts_ws;
        end

        Xstart = [];
        Xstart_norm = [];
        for k = 1:length(sortIdx)
            x = X_cand(sortIdx(k),:);
            x_norm = (x - lb_global) ./ (ub_global - lb_global);
            if isempty(Xstart)
                Xstart = x; Xstart_norm = x_norm;
            else
                d = vecnorm(Xstart_norm - x_norm, 2, 2);
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
            Xstart = [prev_x; Xstart];   % warm start from the previous weight's solution
        end

        objfun = @(x) helper.chebyshevCost(x, Mdl_night, Mdl_day, Mdl_utci, y_min, y_span, lambda, rho);

        best_f = Inf; best_x = []; no_improve = 0;
        for k = 1:size(Xstart,1)
            [x_opt, f_opt] = fmincon(objfun, Xstart(k,:), [],[],[],[], ...
                lb_global, ub_global, [], opts);
            if f_opt < best_f - 1e-8
                best_f = f_opt; best_x = x_opt; no_improve = 0;
            else
                no_improve = no_improve + 1;
            end
            if no_improve >= 5
                break   % stop early once several starts fail to improve
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

%% ===================== FILTER TO NON-DOMINATED =====================
mask = helper.paretoFront(Y_cheby);
Y_cheby_pareto = Y_cheby(mask,:);
X_cheby_pareto = X_cheby(mask,:);

fprintf('%d non-dominated points out of %d Chebyshev samples\n', sum(mask), n_weights);

%% ===================== GP PREDICTIVE UNCERTAINTY ON THE DENSE FRONT =====================
[~, sigmaNight] = predict(Mdl_night, X_cheby_pareto);
[~, sigmaDay]   = predict(Mdl_day,   X_cheby_pareto);
[~, sigmaUTCI]  = predict(Mdl_utci,  X_cheby_pareto);
sigmaTotal = sqrt(sigmaNight.^2 + sigmaDay.^2 + sigmaUTCI.^2);

%% ===================== DOMINANCE CHECK: DO THE LHS PARETO POINTS SURVIVE? =====================
n_true = size(Y_pareto_LHS, 1);
trueStillNonDominated = true(n_true, 1);
for i = 1:n_true
    dominatedBy = all(Y_cheby <= Y_pareto_LHS(i,:), 2) & any(Y_cheby < Y_pareto_LHS(i,:), 2);
    if any(dominatedBy)
        trueStillNonDominated(i) = false;
    end
end

fprintf('%d / %d true Pareto points remain non-dominated after %d GP samples\n', ...
    sum(trueStillNonDominated), n_true, n_weights);

%% ===================== FIGURE: OBJECTIVE-SPACE SCATTER (colored by GP sigma) =====================
fs = 11;

fig = figure('Units','centimeters','Position',[0 0 15 11],'Color','w');
hold on

scatter3(Y_cheby_pareto(:,1), Y_cheby_pareto(:,2), Y_cheby_pareto(:,3), ...
    55, sigmaTotal, 'filled')

scatter3(y_all(:,1), y_all(:,2), y_all(:,3), ...
    20, [0.85 0.7 0.7], 'filled', 'MarkerFaceAlpha', 0.25)

markerSize = 70;
scatter3(Y_pareto_LHS(trueStillNonDominated,1), ...
         Y_pareto_LHS(trueStillNonDominated,2), ...
         Y_pareto_LHS(trueStillNonDominated,3), ...
         markerSize, [0.85 0.1 0.1], 'filled', 'MarkerEdgeColor','k', 'LineWidth',0.6)

scatter3(Y_pareto_LHS(~trueStillNonDominated,1), ...
         Y_pareto_LHS(~trueStillNonDominated,2), ...
         Y_pareto_LHS(~trueStillNonDominated,3), ...
         markerSize, [0.25 0.25 0.25], 'filled', 'MarkerEdgeColor','k', 'LineWidth',0.6)

xlabel('T_{peak}',   'FontName','Times New Roman', 'FontSize',fs, 'Interpreter','tex')
ylabel('T_{night}',  'FontName','Times New Roman', 'FontSize',fs, 'Interpreter','tex')
zlabel('UTCI_{peak}','FontName','Times New Roman', 'FontSize',fs, 'Interpreter','tex')

ax = gca;
ax.FontName  = 'Times New Roman';
ax.FontSize  = fs;
ax.Box       = 'on';
ax.LineWidth = 0.75;
ax.TickDir   = 'out';
grid on
ax.GridAlpha = 0.15;
view(135,25)

colormap(parula)
cb = colorbar;
cb.Label.String   = 'GP predictive \sigma';
cb.Label.FontName = 'Times New Roman';
cb.Label.FontSize = fs;
cb.FontName       = 'Times New Roman';

legend({'GP-sampled front (color = \sigma)', ...
        sprintf('%d LHS samples', m), ...
        sprintf('True non-dominated (still, n=%d)', sum(trueStillNonDominated)), ...
        sprintf('True, now dominated (n=%d)', sum(~trueStillNonDominated))}, ...
    'FontName','Times New Roman', 'FontSize',fs-1, 'Location','northoutside', 'Box','off')

%% ===================== EXPORT FIGURE =====================
if ~exist('figures', 'dir')
    mkdir('figures');
end

exportgraphics(fig, '+figures/lhs_baseline_pareto_densification.png', 'Resolution',600);

%% ===================== SAVE DENSIFICATION RESULTS =====================
% save('+data/+results/LHS_Baseline_Chebyshev.mat', ...
%     'X_cheby', 'Y_cheby', 'X_cheby_pareto', 'Y_cheby_pareto', ...
%     'sigmaNight', 'sigmaDay', 'sigmaUTCI', 'sigmaTotal', ...
%     'trueStillNonDominated', 'n_weights', 'rho', 'minDist');
