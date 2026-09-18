function integrativeslideralpha(DOE,y_night,y_peak,UTCI)
fs = 11;
X = struct2table(DOE);
X.Radius_tree = [];   % drop the derived quantity; alpha stays as the real trained input
set(groot,'defaultAxesFontName','Times')
set(groot,'defaultTextFontName','Times')

tiles   = width(X);
nParams = tiles;

colNames = X.Properties.VariableNames;

% LaTeX display labels — order must match colNames (alpha_tree replaces R_tree)
paramLabels = {'$H_\mathrm{c}$', '$W_\mathrm{c}$', '$W_\mathrm{r}$', ...
               '$f_\mathrm{veg,G}$', '$\alpha_\mathrm{tree}$', '$\alpha_\mathrm{w}$', ...
               '$\lambda_\mathrm{w}$', '$c_\mathrm{v,s,w}$'};

metricLabels = {'$T_\mathrm{peak}$', '$T_\mathrm{night}$', '$\mathrm{UTCI}_\mathrm{peak}$'};
colours = {[0.85 0.2 0.2], [0.2 0.4 0.85], [0.2 0.75 0.35]};

y_night_shape = reshape(y_night,1,[])';
y_peak_shape  = reshape(y_peak,1,[])';
y_UTCI_shape  = reshape(UTCI,1,[])';

rng(16);
n      = height(X);
idx    = randperm(n);
nTrain = round(0.85*n);

trainIdx = idx(1:nTrain);
testIdx  = idx(nTrain+1:end);

Xtrain = X(trainIdx,:);
Xtest  = X(testIdx,:);

y_night_train = y_night_shape(trainIdx);
y_peak_train  = y_peak_shape(trainIdx);
y_utci_train  = y_UTCI_shape(trainIdx);

y_night_test  = y_night_shape(testIdx);
y_peak_test   = y_peak_shape(testIdx);
y_utci_test   = y_UTCI_shape(testIdx);

% =========================================================
% TRAIN GP MODELS (validation split)
% =========================================================
Mdl_day   = fitrgp(Xtrain, y_peak_train,  'KernelFunction','ardsquaredexponential','Standardize',true);
Mdl_night = fitrgp(Xtrain, y_night_train, 'KernelFunction','ardsquaredexponential','Standardize',true);
Mdl_utci  = fitrgp(Xtrain, y_utci_train,  'KernelFunction','ardsquaredexponential','Standardize',true);

% =========================================================
% GP VALIDATION — PARITY PLOTS
% =========================================================
Xtest_tbl = array2table(Xtest{:,:}, 'VariableNames', colNames);

yhat_day   = predict(Mdl_day,   Xtest_tbl);
yhat_night = predict(Mdl_night, Xtest_tbl);
yhat_utci  = predict(Mdl_utci,  Xtest_tbl);

rmse_day   = sqrt(mean((y_peak_test  - yhat_day).^2));
rmse_night = sqrt(mean((y_night_test - yhat_night).^2));
rmse_utci  = sqrt(mean((y_utci_test  - yhat_utci).^2));

r2_day   = 1 - sum((y_peak_test  - yhat_day).^2)   / sum((y_peak_test  - mean(y_peak_test)).^2);
r2_night = 1 - sum((y_night_test - yhat_night).^2) / sum((y_night_test - mean(y_night_test)).^2);
r2_utci  = 1 - sum((y_utci_test  - yhat_utci).^2)  / sum((y_utci_test  - mean(y_utci_test)).^2);

fprintf('RMSE Peak:  %.4f C\n', rmse_day)
fprintf('RMSE Night: %.4f C\n', rmse_night)
fprintf('RMSE UTCI:  %.4f C\n', rmse_utci)

fVal = figure('Name','GP Validation','NumberTitle','off');
set(fVal,'Units','centimeters','Position',[5 5 24 8]);
tVal = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

y_true_all = {y_peak_test,  y_night_test,  y_utci_test};
y_pred_all = {yhat_day,     yhat_night,    yhat_utci};
rmse_all   = {rmse_day,     rmse_night,    rmse_utci};
r2_all     = {r2_day,       r2_night,      r2_utci};

for m = 1:3
    axV = nexttile;
    hold(axV,'on'); box(axV,'on');

    sc = scatter(axV, y_true_all{m}, y_pred_all{m}, 40, colours{m}, 'filled');
    sc.MarkerFaceAlpha = 0.7;

    allVals = [y_true_all{m}; y_pred_all{m}];
    lims    = [min(allVals)-0.5, max(allVals)+0.5];
    plot(axV, lims, lims, 'k--', 'LineWidth', 1);

    xlabel(axV, 'Simulated [°C]',    'FontSize', fs, 'Interpreter','latex');
    ylabel(axV, 'GP predicted [°C]', 'FontSize', fs, 'Interpreter','latex');
    title(axV,  metricLabels{m}, 'FontSize', fs, 'FontWeight','bold', 'Interpreter','latex');
    xlim(axV, lims); ylim(axV, lims);
    axis(axV, 'square');
    set(axV, 'FontSize', fs, 'TickLabelInterpreter','latex');

    text(axV, lims(1)+0.05*diff(lims), lims(2)-0.05*diff(lims), ...
        sprintf('$RMSE = %.3f\\,^{\\circ}C$\n$R^2 = %.3f$', rmse_all{m}, r2_all{m}), ...
        'FontSize', fs, 'VerticalAlignment','top', 'Interpreter','latex');

    grid(axV,'on');
end

exportgraphics(fVal,'GP_validation_alpha.pdf','ContentType','vector');

% =========================================================
% RETRAIN ON FULL DATA — single surrogate, used for slider & Sobol alike
% =========================================================
Mdl_day   = fitrgp(X, y_peak_shape,  'KernelFunction','ardsquaredexponential','Standardize',true);
Mdl_night = fitrgp(X, y_night_shape, 'KernelFunction','ardsquaredexponential','Standardize',true);
Mdl_utci  = fitrgp(X, y_UTCI_shape,  'KernelFunction','ardsquaredexponential','Standardize',true);


% =========================================================
% SOBOL SENSITIVITY — TOTAL ORDER ONLY (Saltelli 2010)
% =========================================================
N     = 16000;
nBoot = 200;

lb = min(X{:,:});
ub = max(X{:,:});

sob        = sobolset(2*nParams, 'Skip',1e3, 'Leap',1e2);
sob        = scramble(sob,'MatousekAffineOwen');
rawSamples = net(sob, N);

A = rawSamples(:, 1:nParams)     .* (ub - lb) + lb;
B = rawSamples(:, nParams+1:end) .* (ub - lb) + lb;

A_tbl = array2table(A, 'VariableNames', colNames);
B_tbl = array2table(B, 'VariableNames', colNames);

yA_day   = predict(Mdl_day,   A_tbl);
yA_night = predict(Mdl_night, A_tbl);
yA_utci  = predict(Mdl_utci,  A_tbl);

yCi_day_all   = zeros(N, nParams);
yCi_night_all = zeros(N, nParams);
yCi_utci_all  = zeros(N, nParams);

for i = 1:nParams
    Ci      = A;
    Ci(:,i) = B(:,i);
    Ci_tbl  = array2table(Ci, 'VariableNames', colNames);
    yCi_day_all(:,i)   = predict(Mdl_day,   Ci_tbl);
    yCi_night_all(:,i) = predict(Mdl_night, Ci_tbl);
    yCi_utci_all(:,i)  = predict(Mdl_utci,  Ci_tbl);
end

ST_day_boot   = zeros(nBoot, nParams);
ST_night_boot = zeros(nBoot, nParams);
ST_utci_boot  = zeros(nBoot, nParams);

for b = 1:nBoot
    bootIdx = randi(N, N, 1);
    for i = 1:nParams
        ST_day_boot(b,i)   = saltelliST(yA_day(bootIdx),   yCi_day_all(bootIdx,i));
        ST_night_boot(b,i) = saltelliST(yA_night(bootIdx), yCi_night_all(bootIdx,i));
        ST_utci_boot(b,i)  = saltelliST(yA_utci(bootIdx),  yCi_utci_all(bootIdx,i));
    end
end

ST_day   = max(mean(ST_day_boot),   0);
ST_night = max(mean(ST_night_boot), 0);
ST_utci  = max(mean(ST_utci_boot),  0);

ST_day_ci   = prctile(ST_day_boot,   [2.5 97.5]);
ST_night_ci = prctile(ST_night_boot, [2.5 97.5]);
ST_utci_ci  = prctile(ST_utci_boot,  [2.5 97.5]);

ST_all    = {ST_day,    ST_night,    ST_utci};
ST_ci_all = {ST_day_ci, ST_night_ci, ST_utci_ci};

% =========================================================
% SENSITIVITY FIGURE — ST ONLY
% =========================================================
fSens = figure('Name','Sobol Sensitivity Analysis','NumberTitle','off');
set(fSens,'Units','centimeters','Position',[5 18 26 10]);

tSens = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

for m = 1:3
    axS = nexttile;
    hold(axS,'on');

    y_pos = (1:nParams)';

    [ST_sorted, sortIdx] = sort(ST_all{m}, 'ascend');
    labels_sorted = paramLabels(sortIdx);

    ST_lo = max(ST_sorted - ST_ci_all{m}(1,sortIdx), 0);
    ST_hi = ST_ci_all{m}(2,sortIdx) - ST_sorted;

    barh(axS, y_pos, ST_sorted, 0.6, ...
        'FaceColor', colours{m}, 'FaceAlpha', 0.75, 'EdgeColor', colours{m});

    errorbar(axS, ST_sorted, y_pos, ST_lo, ST_hi, ...
        'horizontal', 'k.', 'LineWidth', 0.8, 'HandleVisibility','off');

    set(axS, 'YTick', y_pos, 'YTickLabel', labels_sorted, ...
        'FontSize', fs, 'TickLabelInterpreter','latex');
    xlabel(axS, '$S_T$', 'FontSize', fs, 'Interpreter','latex');
    title(axS, metricLabels{m}, 'FontSize', fs, 'FontWeight','bold', 'Interpreter','latex');
    xlim(axS, [0 min(1.05, max(ST_sorted + ST_hi)*1.2)]);
    ylim(axS, [0.5, nParams+0.5]);
    grid(axS,'on');
    box(axS,'on');
end

exportgraphics(fSens,'GP_sobol_ST_alpha.pdf','ContentType','vector');

% =========================================================
% INTERACTIVE SLIDER FIGURE — with live R_tree readout
% =========================================================
f = interactiveGPCutsSlidersAlpha(X, tiles, Mdl_day, Mdl_night, Mdl_utci, ...
    colNames, paramLabels, metricLabels, fs);
drawnow
exportgraphics(f,'GP_slider_snapshot_alpha.pdf','ContentType','vector');

end

% =========================================================
% SALTELLI ST ESTIMATOR
% =========================================================
function ST = saltelliST(yA, yCi)
    % Saltelli et al. (2010) doi:10.1016/j.cpc.2009.09.018
    V = var(yA);
    if V < 1e-12
        ST = 0; return
    end
    ST = mean((yA - yCi).^2) / (2*V);
end




% =========================================================
% INTERACTIVE SLIDER FIGURE (alpha version, with R_tree readout)
% =========================================================
function f = interactiveGPCutsSlidersAlpha(X, tiles, Mdl_day, Mdl_night, Mdl_utci, ...
        colNames, paramLabels, metricLabels, fs)

    fixed = varfun(@mean, X);

    idxAlpha = find(strcmp(colNames, 'alpha'));
    idxWc    = find(strcmp(colNames, 'Width_canyon'));

    f = figure('Name','Interactive GP Cuts with Sliders','NumberTitle','off');
    set(f,'Units','centimeters','Position',[5 5 26 11]);
    tiledlayout(1,tiles,'TileSpacing','compact','Padding','compact');
    ax = gobjects(1,tiles);

    lineDay   = gobjects(1,tiles);
    lineNight = gobjects(1,tiles);
    lineUTCI  = gobjects(1,tiles);
    pCIDay    = gobjects(1,tiles);
    pCINight  = gobjects(1,tiles);
    pCIUTCI   = gobjects(1,tiles);
    vLine     = gobjects(1,tiles);

    for i = 1:tiles
        ax(i) = nexttile;
        xi    = linspace(min(X{:,i}), max(X{:,i}), 100);

        Xquery      = repmat(fixed{1,:}, 100, 1);
        Xquery(:,i) = xi;
        XqueryTbl   = array2table(Xquery, 'VariableNames', colNames);

        [y_day,~,yint_day]     = predict(Mdl_day,   XqueryTbl);
        [y_night,~,yint_night] = predict(Mdl_night, XqueryTbl);
        [y_utci,~,yint_utci]   = predict(Mdl_utci,  XqueryTbl);

        hold(ax(i),'on'); box(ax(i),'on');

        pCIDay(i)   = fill(ax(i),[xi fliplr(xi)], [yint_day(:,1)'   fliplr(yint_day(:,2)')],  'r','FaceAlpha',0.15,'EdgeColor','none');
        pCINight(i) = fill(ax(i),[xi fliplr(xi)], [yint_night(:,1)' fliplr(yint_night(:,2)')],'b','FaceAlpha',0.15,'EdgeColor','none');
        pCIUTCI(i)  = fill(ax(i),[xi fliplr(xi)], [yint_utci(:,1)'  fliplr(yint_utci(:,2)')] ,'g','FaceAlpha',0.15,'EdgeColor','none');

        lineDay(i)   = plot(ax(i), xi, y_day,   'r-',  'LineWidth', 1.5);
        lineNight(i) = plot(ax(i), xi, y_night, 'b--', 'LineWidth', 1.5);
        lineUTCI(i)  = plot(ax(i), xi, y_utci,  'g:',  'LineWidth', 2.0);

        yl   = [12 37];
        xmid = mean([xi(1) xi(end)]);
        vLine(i) = plot(ax(i),[xmid xmid],yl,'k--','LineWidth',1);

        set(ax(i), 'FontSize', fs, 'TickLabelInterpreter','latex');
        title(ax(i), paramLabels{i}, 'Interpreter','latex', 'FontSize', fs)
        xlim(ax(i),[xi(1) xi(end)])
        ylim(ax(i),yl)
        grid(ax(i),'on')

        if strcmp(colNames{i}, 'fveg_G')
            xticks(ax(i), [0.2 0.5 0.8]);
        elseif strcmp(colNames{i}, 'Width_canyon')
            xticks(ax(i), [16 22 28]);
        elseif strcmp(colNames{i}, 'Height_canyon')
            xticks(ax(i), [16 20 24]);
        elseif strcmp(colNames{i}, 'alpha')
            rTreeVal = fixed{1,idxWc} * fixed{1,idxAlpha};
            rTreeAnn = xlabel(ax(i), sprintf('$R_\\mathrm{tree} = %.3f$ m', rTreeVal), ...
                'Interpreter','latex', 'FontSize', fs);
        end

  
    end

    ylabel(ax(1),'Temperature / UTCI [°C]','FontSize',fs,'Interpreter','latex')

    lgd = legend([lineDay(1) lineNight(1) lineUTCI(1) vLine(1)], ...
        [metricLabels, {'Slider position'}], ...
        'Orientation','horizontal', 'FontSize', fs, 'Interpreter','latex');
    lgd.Layout.Tile = 'south';

   

    sliderHeight = 0.02;
    margin       = 0.02;
    usableWidth  = 1 - 2*margin;
    sliderWidth  = usableWidth / tiles;

    for i = 1:tiles
        xpos = margin + (i-1)*sliderWidth;
        uicontrol('Style','slider',...
            'Min',   min(X{:,i}),...
            'Max',   max(X{:,i}),...
            'Value', fixed{1,i},...
            'Units','normalized',...
            'Position',[xpos 0.01 sliderWidth*0.9 sliderHeight],...
            'Callback',@(src,~)updatePlots(src,i));
    end

    function updatePlots(src,dim)
        fixed{1,dim} = src.Value;
        for j = 1:tiles
            xi          = linspace(min(X{:,j}),max(X{:,j}),100);
            Xquery      = repmat(fixed{1,:},100,1);
            Xquery(:,j) = xi;
            XqueryTbl   = array2table(Xquery,'VariableNames',colNames);

            [y_day,~,yint_day]     = predict(Mdl_day,   XqueryTbl);
            [y_night,~,yint_night] = predict(Mdl_night, XqueryTbl);
            [y_utci,~,yint_utci]   = predict(Mdl_utci,  XqueryTbl);

            set(lineDay(j),  'YData', y_day)
            set(lineNight(j),'YData', y_night)
            set(lineUTCI(j), 'YData', y_utci)

            set(pCIDay(j),  'YData', [yint_day(:,1)'   fliplr(yint_day(:,2)')])
            set(pCINight(j),'YData', [yint_night(:,1)' fliplr(yint_night(:,2)')])
            set(pCIUTCI(j), 'YData', [yint_utci(:,1)'  fliplr(yint_utci(:,2)')])

            set(vLine(j),'XData',[fixed{1,j} fixed{1,j}])
      
        end

         % recompute R_tree from current alpha & W_c, update the alpha tile's xlabel
        rTreeVal = fixed{1,idxWc} * fixed{1,idxAlpha};
        set(rTreeAnn, 'String', sprintf('$R_\\mathrm{tree} = %.3f$ m', rTreeVal));
    end

end