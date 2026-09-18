function out = runParEGO_BO(runIdx, DOE0, Params, MeteoData, Name_Site, ...
                             y_night0, y_peak0, y_UTCI_peak0, Results2m0, UTCI0)
%RUNPAREGO_BO Run one independent ParEGO-style Bayesian optimization chain.
%   OUT = RUNPAREGO_BO(RUNIDX, DOE0, PARAMS, METEODATA, NAME_SITE,
%   Y_NIGHT0, Y_PEAK0, Y_UTCI_PEAK0, RESULTS2M0, UTCI0) runs a single
%   independent sequential Bayesian optimization chain, starting from
%   the shared initial LHS archive (DOE0 / y_*0 / RESULTS2M0 / UTCI0),
%   for use as one of the parallel robustness/stability runs. Each
%   iteration: draws a Sobol weight vector on the objective simplex,
%   fits a GP surrogate to the augmented-Tchebycheff (ParEGO)
%   scalarized cost, maximizes Expected Improvement (multi-start
%   fmincon) to propose a new design point, evaluates it with UT&C,
%   and appends it to the archive. Stops on a hard iteration budget or
%   once the hypervolume gain over a trailing window falls below
%   tolerance.
%
%   Inputs:
%       runIdx      - integer run index (1..n_runs); seeds both rng()
%                     and the per-run Sobol stream, so runs are
%                     independent but individually reproducible
%       DOE0        - initial LHS design (struct array), from
%                     generate_initial_LHS_dataset.m
%       Params      - parameter bounds struct (fields .lb / .ub)
%       MeteoData   - meteorological forcing data for UTC.m
%       Name_Site   - site identifier string for UTC.m
%       y_night0, y_peak0, y_UTCI_peak0
%                   - [m x 1] initial LHS objective values
%       Results2m0, UTCI0
%                   - initial LHS simulator output; the archive that
%                     new evaluations get appended to
%
%   Outputs (struct):
%       out.hv_history       - hypervolume at each iteration
%       out.lambda_history    - Sobol simplex weight used each iteration
%       out.y_all              - full archive [T_night, T_peak, UTCI_peak]
%       out.pareto_idx          - indices of out.y_all on the non-dominated front
%       out.DOE, out.X           - full archive design (struct array / table)
%       out.UTCI, out.temp_times - full archive simulator output
%       out.ref_point_final       - fixed reference point used for hypervolume
%       out.n_iters_used           - number of iterations actually run
%       out.final_HV                - hypervolume at the stopping iteration
%       out.pareto_size               - number of non-dominated archive points
%       out.seed                       - rng seed used for this run
%       out.tol_hv, out.hv_window, out.min_iter_before_stop
%                                        - stopping-criterion settings used

ref_point = getFixedRefPoint();   % same reference point used for every run/front

% ---- randomness: independent but reproducible ----
rng(1000 + runIdx);   % distinct, reproducible seed per run

DOE = DOE0; y_night = y_night0; y_peak = y_peak0; y_UTCI_peak = y_UTCI_peak0;
Results2m = Results2m0; UTCI = UTCI0;
y_all = [y_night, y_peak, y_UTCI_peak];

paramNames   = fieldnames(Params);
gpParamNames = paramNames(~strcmp(paramNames, 'Radius_tree'));
nGP = numel(gpParamNames);

X = struct2table(DOE);
if any(strcmp(X.Properties.VariableNames, 'Radius_tree'))
    X.Radius_tree = [];
end

% NOTE: this acquisition-stage rho (0.05) is intentionally different
% from the rho=0.1 used for the post-hoc Pareto-front densification
% scripts (dense_GP_pareto_all_runs.m and
% lhs_baseline_and_chebyshev_densification.m). The two scalarize
% different things - a sequential acquisition criterion here, versus a
% one-shot front reconstruction there - and were tuned separately; this
% is intentional, not a typo.
rho = 0.05;

% ---- distinct Sobol stream per run (different Skip so runs don't reuse the same points) ----
sobolWeights = sobolset(2, 'Skip', 1e3 + runIdx*5000, 'Leap', 1e2);
sobolWeights = scramble(sobolWeights, 'MatousekAffineOwen');

n_iter               = 500;
tol_hv               = 0.001;
hv_window            = 20;
min_iter_before_stop = 50;

hv_history     = [];
lambda_history = [];

lb_global = zeros(1, nGP);
ub_global = zeros(1, nGP);
for j = 1:nGP
    p = gpParamNames{j};
    lb_global(j) = Params.(p).lb;
    ub_global(j) = Params.(p).ub;
end

pareto_idx = [];   % populated on the first iteration; predeclared for scope clarity

for iter = 1:n_iter

    %% ---- Sobol simplex weight for this iteration ----
    uvec = sobolWeights(iter,:);
    u1 = uvec(1); u2 = uvec(2);
    lambda = [1 - sqrt(u1), sqrt(u1)*(1-u2), sqrt(u1)*u2];
    lambda_history = [lambda_history; lambda];

    %% ---- Scalarize the current archive (augmented Tchebycheff) ----
    y_min  = min(y_all, [], 1);
    y_max  = max(y_all, [], 1);
    y_span = max(y_max - y_min, 1e-9);
    y_norm = (y_all - y_min) ./ y_span;

    abs_y_norm    = abs(y_norm);
    weighted_diff = abs_y_norm .* lambda;
    chebyshev_term = max(weighted_diff, [], 2);
    aug_term       = rho * sum(weighted_diff, 2);
    J_cost         = chebyshev_term + aug_term;

    %% ---- Fit GP surrogate of the scalarized cost ----
    GP = fitrgp(X, J_cost, 'KernelFunction','ardsquaredexponential', 'Standardize', true);

    %% ---- Candidate cloud + Expected Improvement ----
    n_candidates = 50000;
    p_sob = sobolset(nGP);
    p_sob = scramble(p_sob, 'MatousekAffineOwen');
    X_cand_unit = net(p_sob, n_candidates);
    X_cand = lb_global + X_cand_unit .* (ub_global - lb_global);

    [mu, sigma] = predict(GP, X_cand);
    sigma = max(sigma, 1e-8);
    f_best = min(J_cost);
    improvement = f_best - mu;
    Z = improvement ./ sigma;
    EI = improvement .* normcdf(Z) + sigma .* normpdf(Z);
    EI(sigma < 1e-8) = 0;

    %% ---- Multi-start selection (greedy, mutually spaced) ----
    nStart = 10;
    minDist = 0.10;
    [~, sortIdx] = sort(EI, 'descend');
    Xstart = []; Xstart_norm = [];
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
        if size(Xstart,1) == nStart
            break
        end
    end

    %% ---- Maximize EI via multi-start fmincon ----
    EIfun = @(x) -helper.expectedImprovementPoint(x, GP, f_best, X.Properties.VariableNames);
    options = optimoptions('fmincon', 'Display','off', 'Algorithm','interior-point', ...
        'StepTolerance',1e-10, 'OptimalityTolerance',1e-8, ...
        'ConstraintTolerance',1e-10, 'MaxFunctionEvaluations',5000);

    bestEI = -Inf; x_new = [];
    duplicateTol = 1e-6;
    for k = 1:size(Xstart,1)
        [xOpt, fval] = fmincon(EIfun, Xstart(k,:), [],[],[],[], ...
            lb_global, ub_global, [], options);
        EIopt = -fval;
        dist = vecnorm(table2array(X) - xOpt, 2, 2);
        if min(dist) < duplicateTol
            continue   % skip points that duplicate an existing archive entry
        end
        if EIopt > bestEI
            bestEI = EIopt; x_new = xOpt;
        end
    end

    if isempty(x_new)
        % fmincon found nothing usable (e.g. every optimum duplicated an
        % existing point); fall back to the best EI candidate directly.
        [~, idx] = max(EI);
        x_new = X_cand(idx,:);
    end

    %% ---- Evaluate the new design point with UT&C ----
    s_new = helper.vec2struct_full(x_new, gpParamNames);
    [result_new, UTCI_new] = UTC(s_new, 1, MeteoData, Name_Site);
    result_new.T2m = result_new.T2m - 273.15;
    Results2m.T2m = cat(3, Results2m.T2m, result_new.T2m);
    UTCI = cat(3, UTCI, UTCI_new);

    y_UTCI_new_peak = max(UTCI_new(:));
    T_peak_new  = max(result_new.T2m(:,:,:), [], 1);
    y_peak_new  = reshape(T_peak_new, 1, [])';
    T_night_new = min(result_new.T2m(40:60,:,:), [], 1);
    y_night_new = reshape(T_night_new, 1, [])';
    y_new = [y_night_new, y_peak_new, y_UTCI_new_peak];

    y_UTCI_peak = [y_UTCI_peak; y_UTCI_new_peak];
    y_peak      = [y_peak; y_peak_new];
    y_night     = [y_night; y_night_new];
    y_all       = [y_all; y_new];

    x_new_table = array2table(x_new, 'VariableNames', X.Properties.VariableNames);
    X = [X; x_new_table];
    DOE(end+1) = s_new;

    %% ---- Update Pareto front / hypervolume, check stopping criterion ----
    isPareto     = helper.paretoFront(y_all);
    pareto_idx   = find(isPareto);
    front_points = y_all(pareto_idx,:);
    HV = helper.hypervolume3D(front_points, ref_point);
    hv_history = [hv_history; HV];

    if iter >= min_iter_before_stop && iter >= hv_window
        hv_prev = hv_history(end-hv_window+1);
        rel_improvement = (hv_history(end) - hv_prev) / max(hv_prev, 1e-9);
        if rel_improvement < tol_hv
            break
        end
    end
end

out.hv_history            = hv_history;
out.lambda_history        = lambda_history;
out.y_all                 = y_all;
out.pareto_idx            = pareto_idx;
out.DOE                   = DOE;
out.X                     = X;
out.UTCI                  = UTCI;
out.temp_times            = Results2m.T2m;
out.ref_point_final       = ref_point;
out.n_iters_used          = numel(hv_history);
out.final_HV              = hv_history(end);
out.pareto_size           = numel(pareto_idx);
out.seed                  = 1000 + runIdx;
out.tol_hv                = tol_hv;
out.hv_window             = hv_window;
out.min_iter_before_stop  = min_iter_before_stop;
end
