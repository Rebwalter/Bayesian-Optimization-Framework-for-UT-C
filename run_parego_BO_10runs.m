%% PARALLEL MULTI-RUN PAREGO BAYESIAN OPTIMIZATION (10 INDEPENDENT RUNS)
%
% Runs 10 independent ParEGO-style Bayesian optimization chains in
% parallel (see runParEGO_BO.m for the per-run algorithm), each
% starting from the same initial LHS archive produced by
% generate_initial_LHS_dataset.m. Used for the stability/robustness
% analysis in the paper: how consistent are the final hypervolume,
% iteration count, and Pareto front size across independent runs.
%
% Requires: +data/+results/Initial_LHS_Dataset.mat (from
%           generate_initial_LHS_dataset.m), runParEGO_BO.m and its
%           dependencies (paretoFront.m, hypervolume3D.m,
%           area2D_dominated.m, expectedImprovementPoint.m,
%           vec2struct_full.m, getFixedRefPoint.m), UTC.m, Parallel
%           Computing Toolbox, Optimization Toolbox, Statistics and
%           Machine Learning Toolbox
%
% Outputs:
%   +data/+results/MOBO_Dataset.mat
%   +data/+results/run_comparison_table.csv
%   figures/hv_convergence_10runs.pdf / .png
%   figures/pareto_fronts_overlay_10runs.pdf / .png

clear; clc;

%% ===================== LOAD INITIAL ARCHIVE =====================
load('+data/+results/Initial_LHS_Dataset.mat')
% Loads DOE, Params, MeteoData, Name_Site, y_all, y_night, y_peak,
% y_UTCI_peak, Results2m, UTCI

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

%% ===================== COMPARE RESULTS ACROSS RUNS =====================
% Single failure-aware summary: runs that hit the catch block above
% have an 'error' field instead of the usual result fields, so every
% aggregate/plot below explicitly excludes them via the `ok` mask
% rather than assuming all runs succeeded.

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

ok = ~failed;

%% ---- Summary table ----
T_compare = table(runIdx, n_iters, final_HV, pareto_size, failed, ...
    'VariableNames', {'Run','Iterations','FinalHV','ParetoSize','Failed'});

disp(T_compare)

%% ---- Aggregate statistics (excluding failed runs) ----
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
[~, rankByHV]    = sort(final_HV(ok), 'descend');
[~, rankByIters] = sort(n_iters(ok));

runsOK = runIdx(ok);
fprintf('\nRuns ranked by final HV (best first): %s\n', mat2str(runsOK(rankByHV)'));
fprintf('Runs ranked by iterations (fewest first): %s\n', mat2str(runsOK(rankByIters)'));

% Run closest to the median final HV (a representative "typical" run)
[~, medianIdxLocal] = min(abs(final_HV(ok) - median(final_HV(ok))));
medianRun = runsOK(medianIdxLocal);
fprintf('Run closest to median final HV: Run %d (HV = %.4f)\n', medianRun, final_HV(medianRun));

%% ---- Export table for appendix ----
if ~exist('+data/+results', 'dir')
    mkdir('+data/+results');
end
writetable(T_compare, '+data/+results/run_comparison_table.csv')

%% ===================== FIGURE: HYPERVOLUME CONVERGENCE OVERLAY =====================
fs = 11;

fig1 = figure('Units','centimeters','Position',[0 0 15 11],'Color','w');
hold on
for r = find(ok)'
    plot(1:numel(results{r}.hv_history), results{r}.hv_history, '-', 'LineWidth', 1)
end
xlabel('BO iteration',  'FontName','Times New Roman', 'FontSize',fs, 'Interpreter','tex')
ylabel('Hypervolume',   'FontName','Times New Roman', 'FontSize',fs, 'Interpreter','tex')

ax = gca;
ax.FontName  = 'Times New Roman';
ax.FontSize  = fs;
ax.Box       = 'on';
ax.LineWidth = 0.75;
ax.TickDir   = 'out';
grid on
ax.GridAlpha = 0.15;

if ~exist('figures', 'dir')
    mkdir('figures');
end
exportgraphics(fig1, '+figures/hv_convergence_10runs.pdf', 'ContentType','vector');
exportgraphics(fig1, '+figures/hv_convergence_10runs.png', 'Resolution',600);

%% ===================== FIGURE: FINAL PARETO FRONTS OVERLAY =====================
colors  = turbo(n_runs);
markers = {'o','s','^','d','v','p','h','o','s','^'};

fig2 = figure('Units','centimeters','Position',[0 0 15 11],'Color','w');
hold on
for r = find(ok)'
    Yp = results{r}.y_all(results{r}.pareto_idx,:);
    scatter3(Yp(:,1), Yp(:,2), Yp(:,3), 35, colors(r,:), markers{r}, 'filled', ...
        'DisplayName', sprintf('Run %d', r))
end
xlabel('T_{night}',  'FontName','Times New Roman', 'FontSize',fs, 'Interpreter','tex')
ylabel('T_{peak}',   'FontName','Times New Roman', 'FontSize',fs, 'Interpreter','tex')
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

legend('show', 'FontName','Times New Roman', 'FontSize',fs-1, 'Location','eastoutside', 'Box','off')

exportgraphics(fig2, '+figures/pareto_fronts_overlay_10runs.png', 'Resolution',600);

%% ===================== SAVE =====================
% save('+data/+results/MOBO_Dataset.mat', 'results', 'n_runs', '-v7.3');
