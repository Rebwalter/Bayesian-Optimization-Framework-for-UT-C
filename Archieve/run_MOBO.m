%% ===================== INITIAL SETUP =====================
load('+data/+results/Initial_LHS_Dataset.mat')
% Loading DOE, Params, MeteoData, Name_Site, y_all, y_night,
% y_peak, y_UTCI_peak, Results2m, UTCI


% ref_point is hard coded to keep consistent


%% ===================== PARALLEL MULTI-RUN DRIVER =====================
if isempty(gcp('nocreate'))
    parpool('local', 10);   % one worker per run; leaves 2 cores free
end

n_runs = 10; 
results = cell(n_runs,1);

parfor r = 1:n_runs
    try
        results{r} = runParEGO_BO(r, DOE, Params, MeteoData, Name_Site, ...
                                   y_night, y_peak, y_UTCI_peak, Results2m, UTCI);
    catch ME
        warning('Run %d failed: %s', r, ME.message);
        results{r} = struct('error', ME.message, 'runIdx', r);
    end
end

%% ===================== COLLECT / COMPARE ACROSS RUNS =====================
final_HV       = cellfun(@(s) s.hv_history(end), results);
n_iters_used   = cellfun(@(s) numel(s.hv_history), results);
pareto_sizes   = cellfun(@(s) numel(s.pareto_idx), results);

fprintf('Final hypervolume across runs: mean=%.4f, CV=%.2f%%\n', ...
    mean(final_HV), 100*std(final_HV)/mean(final_HV));
fprintf('Iterations to stop: mean=%.1f, CV=%.2f%%\n', ...
    mean(n_iters_used), 100*std(n_iters_used)/mean(n_iters_used));
fprintf('Pareto front size: mean=%.1f, CV=%.2f%%\n', ...
    mean(pareto_sizes), 100*std(pareto_sizes)/mean(pareto_sizes));

%% ---- Hypervolume convergence overlay ----
figure; hold on
for r = 1:n_runs
    plot(1:numel(results{r}.hv_history), results{r}.hv_history, '-', 'LineWidth', 1)
end
xlabel('BO Iteration'); ylabel('Hypervolume')
title('Hypervolume convergence across 10 independent BO runs')
grid on

%% ---- Overlay final Pareto fronts across runs ----
figure; hold on
colors = turbo(n_runs);   
markers = {'o','s','^','d','v','p','h','o','s','^'};

for r = 1:n_runs
    Yp = results{r}.y_all(results{r}.pareto_idx,:);
    scatter3(Yp(:,1), Yp(:,2), Yp(:,3), 35, colors(r,:), markers{r}, 'filled', ...
        'DisplayName', sprintf('Run %d', r))
end
xlabel('T_{night}'); ylabel('T_{peak}'); zlabel('UTCI_{peak}')
legend show; grid on; view(135,25)
title('Pareto fronts across independent BO runs')

% save('+data/+results/MOBO_Dataset.mat','results','n_runs','-v7.3');

%% ===================== COMPARE RESULTS ACROSS RUNS =====================

n_runs = numel(results);

runIdx      = (1:n_runs)';
n_iters     = nan(n_runs,1);
final_HV    = nan(n_runs,1);
pareto_size = nan(n_runs,1);
failed      = false(n_runs,1);

for r = 1:n_runs
    if isfield(results{r},'error')
        failed(r) = true;
        continue
    end
    n_iters(r)     = numel(results{r}.hv_history);
    final_HV(r)    = results{r}.hv_history(end);
    pareto_size(r) = numel(results{r}.pareto_idx);
end

%% ---- Summary table ----
T_compare = table(runIdx, n_iters, final_HV, pareto_size, failed, ...
    'VariableNames', {'Run','Iterations','FinalHV','ParetoSize','Failed'});

disp(T_compare)

%% ---- Aggregate statistics (excluding failed runs) ----
ok = ~failed;

fprintf('\n--- Across %d successful runs (of %d total) ---\n', sum(ok), n_runs);

fprintf('Iterations:  mean=%.1f, std=%.2f, CV=%.2f%%, range=[%d, %d]\n', ...
    mean(n_iters(ok)), std(n_iters(ok)), 100*std(n_iters(ok))/mean(n_iters(ok)), ...
    min(n_iters(ok)), max(n_iters(ok)));

fprintf('Final HV:    mean=%.4f, std=%.4f, CV=%.2f%%, range=[%.4f, %.4f]\n', ...
    mean(final_HV(ok)), std(final_HV(ok)), 100*std(final_HV(ok))/mean(final_HV(ok)), ...
    min(final_HV(ok)), max(final_HV(ok)));

fprintf('Pareto size: mean=%.1f, std=%.2f, CV=%.2f%%, range=[%d, %d]\n', ...
    mean(pareto_size(ok)), std(pareto_size(ok)), 100*std(pareto_size(ok))/mean(pareto_size(ok)), ...
    min(pareto_size(ok)), max(pareto_size(ok)));

%% ---- Rank runs (e.g. to identify a "middle" or "best" run) ----
[~, rankByHV]     = sort(final_HV(ok), 'descend');
[~, rankByIters]  = sort(n_iters(ok));

runsOK = runIdx(ok);
fprintf('\nRuns ranked by final HV (best first): %s\n', mat2str(runsOK(rankByHV)'));
fprintf('Runs ranked by iterations (fewest first): %s\n', mat2str(runsOK(rankByIters)'));

% Median-HV run 
[~, medianIdxLocal] = min(abs(final_HV(ok) - median(final_HV(ok))));
medianRun = runsOK(medianIdxLocal);
fprintf('Run closest to median final HV: Run %d (HV = %.4f)\n', medianRun, final_HV(medianRun));

%% ---- Export table for appendix ----
writetable(T_compare, '+data/+results/run_comparison_table.csv')


function out = runParEGO_BO(runIdx, DOE0, Params, MeteoData, Name_Site, ...
                            y_night0, y_peak0, y_UTCI_peak0, Results2m0, UTCI0)
    
    ref_point = [20.5279   35.1913   41.8699];  % hard coded to keep consistent

    % ---- randomness: independent but reproducible ----
    rng(1000 + runIdx);   % distinct, reproducible seed per run

    DOE = DOE0; y_night = y_night0; y_peak = y_peak0; y_UTCI_peak = y_UTCI_peak0;
    Results2m = Results2m0; UTCI = UTCI0;
    y_all = [y_night, y_peak, y_UTCI_peak];

    paramNames = fieldnames(Params);
    gpMask = ~strcmp(paramNames,'Radius_tree');
    gpParamNames = paramNames(gpMask);
    nGP = numel(gpParamNames);

    X = struct2table(DOE);
    if any(strcmp(X.Properties.VariableNames,'Radius_tree'))
        X.Radius_tree = [];
    end

    rho = 0.05;

    % ---- Distinct Sobol stream per run (different Skip so runs don't reuse the same points) ----
    sobolWeights = sobolset(2,'Skip',1e3 + runIdx*5000,'Leap',1e2);
    sobolWeights = scramble(sobolWeights,'MatousekAffineOwen');

    n_iter = 500;
    tol_hv = 0.001; 
    hv_window = 20;
    min_iter_before_stop = 50;

    hv_history = [];
    lambda_history = [];
   
    for j = 1:nGP
        p = gpParamNames{j};
        lb_global(j) = Params.(p).lb;
        ub_global(j) = Params.(p).ub;
    end

    for iter = 1:n_iter
        

        uvec = sobolWeights(iter,:);
        u1 = uvec(1); u2 = uvec(2);
        lambda = [1 - sqrt(u1), sqrt(u1)*(1-u2), sqrt(u1)*u2];
        lambda_history = [lambda_history; lambda];

        y_min = min(y_all,[],1);
        y_max = max(y_all,[],1);
        y_span = max(y_max - y_min, 1e-9);
        y_norm = (y_all - y_min) ./ y_span;

        diff = abs(y_norm);
        weighted_diff = diff .* lambda;
        chebyshev_term = max(weighted_diff, [], 2);
        aug_term = rho * sum(weighted_diff, 2);
        J_cost = chebyshev_term + aug_term;

        GP = fitrgp(X,J_cost,'KernelFunction','ardsquaredexponential','Standardize',true);

        n_candidates = 50000;
        p_sob = sobolset(nGP);
        p_sob = scramble(p_sob,'MatousekAffineOwen');
        X_cand_unit = net(p_sob,n_candidates);
        X_cand = lb_global + X_cand_unit .* (ub_global-lb_global);

        [mu,sigma] = predict(GP,X_cand);
        sigma = max(sigma,1e-8);
        f_best = min(J_cost);
        improvement = f_best - mu;
        Z = improvement ./ sigma;
        EI = improvement .* normcdf(Z) + sigma .* normpdf(Z);
        EI(sigma < 1e-8) = 0;

        nStart = 10;
        minDist = 0.10;
        [~,sortIdx] = sort(EI,'descend');
        Xstart = []; Xstart_norm = [];
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
            if size(Xstart,1)==nStart, break; end
        end

        EIfun = @(x) -expectedImprovementPoint(x,GP,f_best,X.Properties.VariableNames);
        options = optimoptions('fmincon','Display','off','Algorithm','interior-point', ...
            'StepTolerance',1e-10,'OptimalityTolerance',1e-8, ...
            'ConstraintTolerance',1e-10,'MaxFunctionEvaluations',5000);

        bestEI = -Inf; x_new = [];
        duplicateTol = 1e-6;
        for k = 1:size(Xstart,1)
            [xOpt,fval] = fmincon(EIfun,Xstart(k,:),[],[],[],[], ...
                lb_global,ub_global,[],options);
            EIopt = -fval;
            dist = vecnorm(table2array(X)-xOpt,2,2);
            if min(dist) < duplicateTol, continue; end
            if EIopt > bestEI
                bestEI = EIopt; x_new = xOpt;
            end
        end

        if isempty(x_new)
            [~,idx] = max(EI);
            x_new = X_cand(idx,:);
        end

        s_new = vec2struct_full(x_new,gpParamNames);
        [result_new,UTCI_new] = UTC(s_new,1,MeteoData,Name_Site);
        result_new.T2m = result_new.T2m-273.15;
        Results2m.T2m = cat(3,Results2m.T2m,result_new.T2m);
        UTCI = cat(3,UTCI,UTCI_new);

        y_UTCI_new_peak = max(UTCI_new(:));
        T_peak_new = max(result_new.T2m(:,:,:),[],1);
        y_peak_new = reshape(T_peak_new,1,[])';
        T_night_new = min(result_new.T2m(40:60,:,:),[],1);
        y_night_new = reshape(T_night_new,1,[])';
        y_new = [y_night_new, y_peak_new, y_UTCI_new_peak];

        y_UTCI_peak = [y_UTCI_peak; y_UTCI_new_peak];
        y_peak = [y_peak; y_peak_new];
        y_night = [y_night; y_night_new];
        y_all = [y_all; y_new];

        x_new_table = array2table(x_new,'VariableNames',X.Properties.VariableNames);
        X = [X; x_new_table];
        DOE(end+1) = s_new;

        n_all = size(y_all,1);
        isPareto = true(n_all,1);
        for i = 1:n_all
            for jj = 1:n_all
                if i ~= jj
                    if all(y_all(jj,:) <= y_all(i,:)) && any(y_all(jj,:) < y_all(i,:))
                        isPareto(i) = false;
                        break
                    end
                end
            end
        end
        pareto_idx = find(isPareto);
        front_points = y_all(pareto_idx,:);
        HV = hypervolume3D(front_points, ref_point);
        hv_history = [hv_history; HV];

        if iter >= min_iter_before_stop && iter >= hv_window
            hv_prev = hv_history(end-hv_window+1);
            rel_improvement = (hv_history(end) - hv_prev) / max(hv_prev, 1e-9);
            if rel_improvement < tol_hv
                break
            end
        end
    end

    out.hv_history      = hv_history;
    out.lambda_history  = lambda_history;
    out.y_all           = y_all;
    out.pareto_idx       = pareto_idx;
    out.DOE              = DOE;
    out.X                = X;
    out.UTCI             = UTCI;
    out.temp_times       = Results2m.T2m;
    out.ref_point_final  = ref_point;               
    out.n_iters_used     = numel(hv_history);
    out.final_HV         = hv_history(end);
    out.pareto_size      = numel(pareto_idx);
    out.seed             = 1000 + runIdx;
    out.tol_hv           = tol_hv;
    out.hv_window        = hv_window;
    out.min_iter_before_stop = min_iter_before_stop;
end




%%%% Helper
function EI = expectedImprovementPoint(x,GP,f_best,varNames)

% x is a 1xd vector


[mu,sigma] = predict(GP,x);

sigma = max(sigma,1e-6);

improvement = f_best - mu;
Z = improvement ./ sigma;

EI = improvement .* normcdf(Z) + sigma .* normpdf(Z);

if sigma < 1e-8
    EI = 0;
end

end


function s = vec2struct_full(x,gpParamNames)
    for j = 1:numel(gpParamNames)
        s.(gpParamNames{j}) = x(j);
    end
    s.Radius_tree = s.alpha * s.Width_canyon;
end

function HV = hypervolume3D(P, r)
% 3D hypervolume via dimension-sweep (HSO), minimization convention.
% P: n x 3 non-dominated points. r: 1x3 reference point (worse than all P).
%
% Convention: sort ASCENDING by objective 1 (best/smallest first). The
% slab for point k covers x in [P(k,1), x_next), and is "covered" by
% every point with x1 <= P(k,1) -- i.e. points 1..k in this ordering --
% since each point's dominated box extends from its own x1 out to r(1).
 
    if isempty(P)
        HV = 0;
        return
    end
 
    P = sortrows(P, 1);   % ascending by objective 1
    n = size(P,1);
    HV = 0;
 
    for k = 1:n
        if k < n
            x_next = P(k+1,1);
        else
            x_next = r(1);
        end
        width = x_next - P(k,1);
        if width <= 0
            continue
        end
 
        active = P(1:k, 2:3);
        area = area2D_dominated(active, r(2), r(3));
 
        HV = HV + width * area;
    end
end

function area = area2D_dominated(points, r2, r3)
% 2D dominated area union (skyline sweep), minimization convention.
% Same logic as hypervolume3D one dimension down: sort ascending by y,
% accumulate the running best (smallest) z, and each point's slab runs
% from its own y out to the next point's y (or r2 for the last one).
    points = sortrows(points, 1);   % ascending by y
    n = size(points,1);
 
    area = 0;
    z_min = Inf;
 
    for k = 1:n
        z_min = min(z_min, points(k,2));
 
        if k < n
            y_next = points(k+1,1);
        else
            y_next = r2;
        end
 
        width = y_next - points(k,1);
        height = r3 - z_min;
 
        area = area + width * height;
    end
end
 