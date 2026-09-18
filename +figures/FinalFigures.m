%% ===================== EXTRACT RUN 8 =====================
load('2_BOxTchebysheffx10Parallel.mat')  
load('3_GPbasedParetoFrontsx10.mat') 

result8 = results{8,1};

y_all          = result8.y_all;
pareto_idx     = result8.pareto_idx;
DOE            = result8.DOE;
X              = result8.X;
hv_history     = result8.hv_history;
lambda_history = result8.lambda_history;
UTCI           = result8.UTCI;
Results2m.T2m  = result8.temp_times;

m = 50; 
iter = size(hv_history,1);

DOE_pareto = DOE(pareto_idx);

y_night   = y_all(:,1);
y_peak    = y_all(:,2);
UTCI_peak = y_all(:,3);

Tpeak_min = min(y_all(:,2));
Tpeak_max = max(y_all(:,2));

fprintf('Run 8: %d iterations, final HV = %.4f, %d Pareto points\n', ...
    numel(hv_history), hv_history(end), numel(pareto_idx));

%% ===================== GP DENSE PARETO FRONT (RUN 8) =====================
% Assumes GP_dense_results(runIdx) internally pulls from results{runIdx}
% and returns/assigns Y_cheby_pareto, sigmaTotal into the base workspace.

gpDense8 = GP_dense_results(8);

Y_cheby_pareto = gpDense8.Y_cheby_pareto;
X_cheby_pareto = gpDense8.X_cheby_pareto;
sigmaTotal     = gpDense8.sigmaTotal;
n_pareto_dense = gpDense8.n_pareto;
HV_dense       = gpDense8.HV_dense_fixedref;
HV_archive     = gpDense8.HV_archive_fixedref;
%% ===================== 3-VIEW PARETO SCATTER (RUN 8) =====================
figure('Position',[100 100 1400 420]);
t = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

viewAngles = {[0 90], [0 0], [90 0]};   % top, front (side), right (side)

for p = 1:3
    ax = nexttile; hold on

    scatter3(y_all(:,1), y_all(:,2), y_all(:,3),...
        15,[0.8 0.8 0.8],'filled', 'DisplayName','True dominated (UT&C)')
    scatter3(y_all(pareto_idx,1),...
             y_all(pareto_idx,2),...
             y_all(pareto_idx,3),...
             90,'k','filled',...
             'LineWidth',1.2,...
             'DisplayName','True non-dominated (UT&C)')
    scatter3(Y_cheby_pareto(:,1),...
             Y_cheby_pareto(:,2),...
             Y_cheby_pareto(:,3),...
             70,...
             sigmaTotal,...
             'filled', 'DisplayName','GP sampled Pareto front')

    xlabel('T_{night}')
    ylabel('T_{peak}')
    zlabel('UTCI_{peak}')
    zlim([26.5,35]);
    grid on
    box(ax,'on')
    axis(ax,'square')
    view(viewAngles{p})
    colormap(parula)
end

cb = colorbar;
cb.Layout.Tile = 'east';
cb.Label.String = 'GP predictive \sigma';
clim([0 0.6])

lgd = legend('show','Orientation','horizontal');
lgd.Layout.Tile = 'south';
lgd.Box = 'off';

set(findall(gcf,'-property','FontName'),'FontName','Times New Roman')
set(findall(gcf,'-property','FontSize'),'FontSize',11)
% exportgraphics(gcf,'run8_pareto_3views.pdf','ContentType','vector')
exportgraphics(gcf,'run8_pareto_3views.png','Resolution',600)


%% ===================== PARETO RANGE / PARAMETER DISTRIBUTION (RUN 8) =====================
varNames    = {'height_c','width_c','width_r','fveg_g','alpha','albedo_w','lan_w','cv_w'};
paramLabels = { ...
    '$H_\mathrm{c}$', ...
    '$W_\mathrm{c}$', ...
    '$W_\mathrm{r}$', ...
    '$f_\mathrm{veg,G}$', ...
    '$\alpha_\mathrm{tree}$', ...
    '$\alpha_\mathrm{w}$', ...
    '$\lambda_\mathrm{w}$', ...
    '$c_\mathrm{v,s,w}$'};
paramNames  = {'Height_canyon','Width_canyon','Width_roof','fveg_G','alpha','albedo_w','lan_dry_W','cv_s_W'};

% Highlighted runs (best per-objective points within run 8's dataset)
[~,idxnight] = min(y_all(:,1));
[~,idxpeak]  = min(y_all(:,2));
[~,idxUTCI]  = min(y_all(:,3));

highlightRuns = [idxpeak,idxnight,idxUTCI];

highlightColors = [ ...
    0.85 0.20 0.20;   % best T_peak     (red)
    0.20 0.45 0.90;   % best T_night    (blue)
    0.20 0.70 0.25];  % best utci       (green)
highlightMarkers = {'o','s','^'};

nVar  = length(varNames);
nCols = 8;
nRows = ceil(nVar/nCols);

f = figure;
set(f,'Units','normalized')

t = tiledlayout(nRows, nCols, ...
    'TileSpacing','compact', ...
    'Padding','compact');

lgdHandles = gobjects(6,1);

for i = 1:nVar

    nexttile
    hold on

    varName   = varNames{i};
    paramName = paramNames{i};

    x_par = [DOE_pareto.(paramName)]';

    if strcmp(varName,'Radius_tree')
        x_all = [DOE.(paramName)]';
        lb = min(x_all);
        ub = max(x_all);
    else
        lb = Params.(paramName).lb;
        ub = Params.(paramName).ub;
    end

    if strcmp(varName, 'cv_w')
        x_par = x_par / 1e6;
        lb    = lb    / 1e6;
        ub    = ub    / 1e6;
    end

    hBounds = plot([1 1], [lb ub], ...
        'Color',[0.75 0.75 0.75], ...
        'LineWidth',8);

    plot(1, lb, '_k', 'MarkerSize',18, 'LineWidth',2)
    plot(1, ub, '_k', 'MarkerSize',18, 'LineWidth',2)

    hBox = boxchart(ones(size(x_par)), x_par, ...
        'BoxFaceColor',[0.8 0.2 0.2], ...
        'BoxWidth',0.35, ...
        'MarkerStyle','none');

    parMin = min(x_par);
    parMax = max(x_par);
    cRange = [0.5 0.0 0.0];

    hParetoRange = plot(1, parMin, '_', ...
        'Color',cRange, 'MarkerSize',18, 'LineWidth',3);
    plot(1, parMax, '_', ...
        'Color',cRange, 'MarkerSize',18, 'LineWidth',3)

    for j = 1:length(highlightRuns)

        runID = highlightRuns(j);
        val   = DOE(runID).(paramName);

        if strcmp(varName,'cv_w')
            val = val/1e6;
        end

        h = scatter(1, val, 90, ...
            'filled', ...
            'Marker', highlightMarkers{j}, ...
            'MarkerFaceColor', highlightColors(j,:), ...
            'MarkerEdgeColor', 'k');

        if i == 1
            lgdHandles(j) = h;
        end
    end

    if i == 1
        lgdHandles(4) = hBox;
        lgdHandles(5) = hBounds;
        lgdHandles(6) = hParetoRange;
    end

    title(paramLabels{i}, 'Interpreter','latex')
    set(gca, 'XTick',[], 'FontSize',12)
    xlim([0.7 1.3])

    if strcmp(varName, 'cv_w')
        ylabel('$[\times 10^6]$', 'Interpreter','latex', 'FontSize',9)
    end

    if strcmp(varName, 'fveg_g')
        ylim([0.0 1.05])
    end

    grid on
    box on

end

lgdLabels = { ...
    'Best $T_\mathrm{peak}$', ...
    'Best $T_\mathrm{night}$', ...
    'Best $\mathrm{UTCI}_\mathrm{peak}$', ...
    'Distribution non-dominated', ...
    '  Parameter range', ...
    '  Range non-dominated'};

lgd = legend(lgdHandles, lgdLabels, ...
    'Orientation','horizontal', ...
    'Box','off', ...
    'Interpreter','latex', 'FontSize', 12);
lgd.Layout.Tile = 'south';
lgd.ItemTokenSize = [45 18];

set(findall(f,'-property','FontName'),'FontName','Times New Roman')
exportgraphics(f, 'run8_ParetoRanges.pdf', 'ContentType','vector');





%% ===================== SLIDER (RUN 8) =====================
helper.integrativeslideralpha(DOE,y_night,y_peak,UTCI_peak)

helper.ReadPostprocessing(...
    MeteoData,...
    Results2m.T2m,...
    Results2m.q2m,UTCI,...
    m,iter,...
    MeteoData.Tatm);