function [] = ReadPostprocessing(Meteodata, T2m, q2m, UTCI, m, iter, T_atm)

fs    = 14;
Date  = Meteodata.Time;
nRuns = m + iter;

% Peak temperature per run for color coding
peakT   = squeeze(max(T2m, [], [1 2]));
peakT_n = (peakT - min(peakT)) / (max(peakT) - min(peakT));

% High contrast light pink to deep dark red
cmap_red = [linspace(1.00, 0.50, 256)', ...
            linspace(0.85, 0.00, 256)', ...
            linspace(0.85, 0.00, 256)'];

colIdx     = max(1, round(peakT_n * 255) + 1);
runColours = cmap_red(colIdx, :);

% =========================================================
% TEMPERATURE + UTCI — side by side
% =========================================================
figTU = figure('Name','Temperature and UTCI time series','NumberTitle','off');
set(figTU,'Units','centimeters','Position',[5 7 30 10]);
tTU = tiledlayout(figTU, 1, 2, 'TileSpacing','compact','Padding','compact');

% --- Temperature ---
axT = nexttile(tTU);
hold(axT,'on'); box(axT,'on');

for i = 1:nRuns
    h = plot(axT, Date, T2m(:,:,i), ...
        'Color', [runColours(i,:) 0.5], ...
        'LineWidth', 1.5, ...
        'HandleVisibility','off');
    h.DataTipTemplate.DataTipRows(1).Label = 'Time';
    h.DataTipTemplate.DataTipRows(2).Label = 'Temperature [°C]';
    h.DataTipTemplate.DataTipRows(end+1)   = ...
        dataTipTextRow('Run index', repmat(i, size(Date)));
end

hAtm = plot(axT, Date, T_atm - 273.15, ...
    'Color', [0.15 0.15 0.15], ...
    'LineWidth', 1.7, ...
    'DisplayName', 'T_{atm}');
legend(axT, hAtm, ...
    'Interpreter','tex', ...
    'FontSize', fs, 'FontName','Times New Roman', ...
    'Location','northwest', 'Box','off');
xlabel(axT, 'Date',             'FontSize', fs, 'Interpreter','latex', 'FontName','Times New Roman');
ylabel(axT, 'Temperature [°C]','FontSize', fs, 'Interpreter','latex', 'FontName','Times New Roman');
set(axT, 'FontSize', fs, 'FontName','Times New Roman', ...
    'TickLabelInterpreter','latex', 'LineWidth', 0.8);
xtickformat('MM-dd');
ylim(axT, [12 32]);
grid(axT,'on');
colormap(axT, cmap_red);
clim(axT, [min(peakT) max(peakT)]);

% --- UTCI ---
axU = nexttile(tTU);
hold(axU,'on'); box(axU,'on');

for i = 1:nRuns
    h = plot(axU, Date, UTCI(:,:,i), ...
        'Color', [runColours(i,:) 0.5], ...
        'LineWidth', 1.5, ...
        'HandleVisibility','off');
    h.DataTipTemplate.DataTipRows(1).Label = 'Time';
    h.DataTipTemplate.DataTipRows(2).Label = 'UTCI [°C]';
    h.DataTipTemplate.DataTipRows(end+1)   = ...
        dataTipTextRow('Run index', repmat(i, size(Date)));
end

xlabel(axU, 'Date',      'FontSize', fs, 'Interpreter','latex', 'FontName','Times New Roman');
ylabel(axU, 'UTCI [°C]','FontSize', fs, 'Interpreter','latex', 'FontName','Times New Roman');
set(axU, 'FontSize', fs, 'FontName','Times New Roman', ...
    'TickLabelInterpreter','latex', 'LineWidth', 0.8);
xtickformat('MM-dd');
ylim(axU, [12 37]);
grid(axU,'on');
colormap(axU, cmap_red);
clim(axU, [min(peakT) max(peakT)]);

% --- Shared colorbar on the right ---
cb = colorbar(axU);
clim(axU, [min(peakT) max(peakT)]);
cb.Label.String      = 'T_{peak} [°C]';
cb.Label.Interpreter = 'tex';
cb.Label.FontSize    = fs;
cb.Label.FontName    = 'Times New Roman';
cb.TickLabelInterpreter = 'tex';
cb.FontSize          = fs;

% --- Shared legend between the two panels ---
% lgd = legend(axT, hAtm, ...
%     'Interpreter','tex', ...
%     'FontSize', fs, 'FontName','Times New Roman', ...
%     'Box','off');
% lgd.Layout.Tile = 'south';

exportgraphics(figTU, sprintf('TempUTCI_timeseries_%dRuns.pdf', nRuns), ...
    'ContentType','vector');

% =========================================================
% HUMIDITY
% =========================================================
figQ = figure('Name','Humidity time series','NumberTitle','off');
set(figQ,'Units','centimeters','Position',[5 7 16 10]);
axQ  = axes(figQ);
hold(axQ,'on'); box(axQ,'on');

for i = 1:m
    h = plot(axQ, Date, q2m(:,:,i), ...
        'Color', [runColours(i,:) 0.5], ...
        'LineWidth', 1.5, ...
        'HandleVisibility','off');
    h.DataTipTemplate.DataTipRows(1).Label = 'Time';
    h.DataTipTemplate.DataTipRows(2).Label = 'Humidity [kg/kg]';
    h.DataTipTemplate.DataTipRows(end+1)   = ...
        dataTipTextRow('Run index', repmat(i, size(Date)));
end

colormap(axQ, cmap_red);
cb3 = colorbar(axQ);
clim(axQ, [min(peakT) max(peakT)]);
cb3.Label.String      = 'T_{peak} [°C]';
cb3.Label.Interpreter = 'tex';
cb3.Label.FontSize    = fs;
cb3.Label.FontName    = 'Times New Roman';
cb3.TickLabelInterpreter = 'tex';
cb3.FontSize          = fs;

xlabel(axQ, 'Date',             'FontSize', fs, 'Interpreter','latex', 'FontName','Times New Roman');
ylabel(axQ, 'Humidity [kg/kg]','FontSize', fs, 'Interpreter','latex', 'FontName','Times New Roman');
set(axQ, 'FontSize', fs, 'FontName','Times New Roman', ...
    'TickLabelInterpreter','latex', 'LineWidth', 0.8);
xtickformat('MM-dd');
grid(axQ,'on');

exportgraphics(figQ, sprintf('Humidity_timeseries_%dRuns.pdf', nRuns), ...
    'ContentType','vector');

end