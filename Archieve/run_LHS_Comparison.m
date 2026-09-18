%%%%%%%%%% UTCI MIN / MAX OPTIMISATION (PARALLELIZED) %%%%%
clear; clc;

%% ===================== SETUP =====================
load('MeteoData_DWD_hourly_RadPart.mat');
MeteoData = MeteoData_DWD(1009:1105,:);
Name_Site = 'HST_BEM';
m = 144;
rdmseed = 11; 

%% ===================== INITIAL DOE (unchanged, sequential) =====================
[DOE, ~, Params] = makeDOE_opt(m,rdmseed);

%% ===================== START PARALLEL POOL =====================
if isempty(gcp('nocreate'))
    parpool('local', 11);
end
poolObj  = gcp();
nWorkers = poolObj.NumWorkers;
fprintf('Using %d parallel workers\n', nWorkers);

%% ===================== CHOP DOE INTO CHUNKS =====================
chunkEdges = round(linspace(0, m, nWorkers+1));

DOE_chunks = cell(nWorkers,1);
m_chunks   = zeros(nWorkers,1);

for c = 1:nWorkers
    idxRange = (chunkEdges(c)+1):chunkEdges(c+1);
    DOE_chunks{c} = DOE(idxRange);
    m_chunks(c)   = numel(idxRange);
end

%% ===================== PARALLEL UTC EVALUATION =====================
T2m_chunks  = cell(nWorkers,1);
UTCI_chunks = cell(nWorkers,1);

parfor c = 1:nWorkers

    DOE_c = DOE_chunks{c};
    m_c   = m_chunks(c);

    [Results2m_c, UTCI_c] = UTC(DOE_c, m_c, MeteoData, Name_Site);
    Results2m_c.T2m = Results2m_c.T2m - 273.15;

    T2m_chunks{c}  = Results2m_c.T2m;
    UTCI_chunks{c} = UTCI_c;

end

%% ===================== REASSEMBLE IN ORIGINAL DOE ORDER =====================
Results2m.T2m = cat(3, T2m_chunks{:});
UTCI          = cat(3, UTCI_chunks{:});

%% ===================== SANITY CHECK: reassembled size matches m =====================
assert(size(Results2m.T2m,3) == m, ...
    'Reassembled T2m has %d slices, expected m=%d', size(Results2m.T2m,3), m);
assert(size(UTCI,3) == m, ...
    'Reassembled UTCI has %d slices, expected m=%d', size(UTCI,3), m);

%% ===================== OBJECTIVES (unchanged) =====================
% Daytime peak
T_peak = max(Results2m.T2m(:,:,:), [], 1);
y_peak = reshape(T_peak,1,[])';

% Nighttime low
T_window = Results2m.T2m(40:60,:,:);
NightMinTemp = min(T_window, [], 1);
y_night = reshape(NightMinTemp,1,[])';

% UTCI peak
UTCI_peak = max(UTCI,[],1);
UTCI_peak = reshape(UTCI_peak,1,[])';

y_all = [y_night, y_peak, UTCI_peak];
Y     = [y_night, y_peak, UTCI_peak];

%% ===================== NON-DOMINATED (PARETO) SUBSET =====================
n_all = size(y_all,1);
isPareto = true(n_all,1);
for i = 1:n_all
    for j = 1:n_all
        if i ~= j
            if all(y_all(j,:) <= y_all(i,:)) && any(y_all(j,:) < y_all(i,:))
                isPareto(i) = false;
                break
            end
        end
    end
end
pareto_idx_LHS = find(isPareto);
Y_pareto_LHS   = y_all(pareto_idx_LHS,:);

fprintf('LHS baseline: %d design points, %d non-dominated points\n', ...
    m, numel(pareto_idx_LHS));

%% ===================== HYPERVOLUME (fixed ref point) =====================
ref_point_fixed = [20.5279   35.1913   41.8699];   % same fixed box used everywhere else
HV_LHS = hypervolume3D(Y_pareto_LHS, ref_point_fixed);

fprintf('LHS baseline hypervolume (fixed ref): %.4f\n', HV_LHS);

%% ===================== SAVE =====================
% save('+data\+results\CompareLHS144_random.mat','y_all','Y_pareto_LHS','pareto_idx_LHS','HV_LHS', ...
%     'ref_point_fixed','Name_Site',"Params",'DOE',"m","Results2m",'UTCI', ...
%     'y_peak','y_night','UTCI_peak','MeteoData')
% 
%% ===================== TRAIN 3 GP SURROGATES =====================
paramNames = fieldnames(Params);
gpMask = ~strcmp(paramNames,'Radius_tree');
gpParamNames = paramNames(gpMask);

X = struct2table(DOE);
if any(strcmp(X.Properties.VariableNames,'Radius_tree'))
    X.Radius_tree = [];
end

y_night_shape = y_all(:,1);
y_peak_shape  = y_all(:,2);
y_UTCI_shape  = y_all(:,3);

Mdl_night = fitrgp(X, y_night_shape, 'KernelFunction','ardsquaredexponential','Standardize',true);
Mdl_day   = fitrgp(X, y_peak_shape,  'KernelFunction','ardsquaredexponential','Standardize',true);
Mdl_utci  = fitrgp(X, y_UTCI_shape,  'KernelFunction','ardsquaredexponential','Standardize',true);

fprintf('GPs trained: Mdl_night, Mdl_day, Mdl_utci (on %d points, %d parameters)\n', ...
    size(X,1), size(X,2));

% save('GPs_144_fixedref.mat','Mdl_night','Mdl_day','Mdl_utci','y_pareto','pareto_idx', ...
%     'DOE_pareto','HV_144','ref_point_fixed','X')

%% ===================== HELPER FUNCTIONS =====================
function mask = paretoFront(Y)
    n = size(Y,1);
    mask = true(n,1);
    for i = 1:n
        if mask(i)
            dominates = all(Y <= Y(i,:), 2) & any(Y < Y(i,:), 2);
            dominates(i) = false;
            if any(dominates)
                mask(i) = false;
            end
        end
    end
end

function HV = hypervolume3D(P, r)
    if isempty(P)
        HV = 0;
        return
    end
    P = sortrows(P, 1);
    n = size(P,1);
    HV = 0;
    for k = 1:n
        if k < n
            x_next = P(k+1,1);
        else
            x_next = r(1);
        end
        width = x_next - P(k,1);
        if width <= 0, continue; end
        active = P(1:k, 2:3);
        area = area2D_dominated(active, r(2), r(3));
        HV = HV + width * area;
    end
end

function area = area2D_dominated(points, r2, r3)
    points = sortrows(points, 1);
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


%% ===================== SETUP =====================


if isempty(gcp('nocreate'))
    parpool('local');
end
poolObj = gcp();
nWorkers = poolObj.NumWorkers-1;
fprintf('Using %d parallel workers\n', nWorkers);

paramNames = fieldnames(Params);
gpMask = ~strcmp(paramNames,'Radius_tree');
gpParamNames = paramNames(gpMask);
nGP = numel(gpParamNames);

lb_global = zeros(1,nGP);
ub_global = zeros(1,nGP);
for j = 1:nGP
    p = gpParamNames{j};
    lb_global(j) = Params.(p).lb;
    ub_global(j) = Params.(p).ub;
end

%% ===================== SAMPLE 10000 POINTS VIA CHEBYSHEV + GPs =====================
n_weights   = 10000;
n_starts    = 20;
n_starts_ws = 8;
minDist     = 0.10;
rho         = 0.1;

opts = optimoptions('fmincon', ...
    'Display','off','Algorithm','sqp', ...
    'MaxFunctionEvaluations',500, ...
    'OptimalityTolerance',1e-8, ...
    'StepTolerance',1e-10, ...
    'ConstraintTolerance',1e-8);

X_cand_unit = lhsdesign(10000,nGP);
X_cand = lb_global + (ub_global-lb_global).*X_cand_unit;

y_min  = min(y_all,[],1);
y_max  = max(y_all,[],1);
y_span = max(y_max-y_min,1e-9);

yNight = predict(Mdl_night,X_cand);
yDay   = predict(Mdl_day,X_cand);
yUTCI  = predict(Mdl_utci,X_cand);
Ynorm  = ([yNight yDay yUTCI]-y_min)./y_span;

p = sobolset(2);
p = scramble(p,'MatousekAffineOwen');
sobolWeights = net(p,n_weights);

[sobolWeights_sorted, sortWeightIdx] = sortrows(sobolWeights, [1 2]);
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

        objfun = @(x) chebyshevCost(x, Mdl_night, Mdl_day, Mdl_utci, y_min, y_span, lambda, rho);

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

%% ---- Filter GP samples to non-dominated ----
mask = paretoFront(Y_cheby);
Y_cheby_pareto = Y_cheby(mask,:);
X_cheby_pareto = X_cheby(mask,:);

fprintf('%d non-dominated points out of %d Chebyshev samples\n', sum(mask), n_weights);

%% ---- Sigma on the GP-dense Pareto points ----
[~,sigmaNight] = predict(Mdl_night, X_cheby_pareto);
[~,sigmaDay]   = predict(Mdl_day,   X_cheby_pareto);
[~,sigmaUTCI]  = predict(Mdl_utci,  X_cheby_pareto);
sigmaTotal = sqrt(sigmaNight.^2 + sigmaDay.^2 + sigmaUTCI.^2);

%% ===================== CHECK DOMINANCE OF THE 22 TRUE PARETO POINTS =====================
n_true = size(y_pareto,1);
trueStillNonDominated = true(n_true,1);

for i = 1:n_true
    dominatedBy = all(Y_cheby <= y_pareto(i,:), 2) & any(Y_cheby < y_pareto(i,:), 2);
    if any(dominatedBy)
        trueStillNonDominated(i) = false;
    end
end

fprintf('%d / %d true Pareto points remain non-dominated after 2000 GP samples\n', ...
    sum(trueStillNonDominated), n_true);

%% ===================== PLOT =====================
fs = 11;

figure('Units','centimeters','Position',[0 0 15 11],'Color','w')
hold on

scatter3(Y_cheby_pareto(:,1), Y_cheby_pareto(:,2), Y_cheby_pareto(:,3), ...
    55, sigmaTotal, 'filled')

markerSize = 70;
scatter3(y_all(:,1), ...
         y_all(:,2), ...
         y_all(:,3), ...
         20,[0.85 0.7 0.7],...
         'filled','MarkerFaceAlpha',0.25)

scatter3(y_pareto(trueStillNonDominated,1), ...
         y_pareto(trueStillNonDominated,2), ...
         y_pareto(trueStillNonDominated,3), ...
         markerSize, [0.85 0.1 0.1], 'filled', 'MarkerEdgeColor','k', 'LineWidth',0.6)

scatter3(y_pareto(~trueStillNonDominated,1), ...
         y_pareto(~trueStillNonDominated,2), ...
         y_pareto(~trueStillNonDominated,3), ...
         markerSize, [0.25 0.25 0.25], 'filled', 'MarkerEdgeColor','k', 'LineWidth',0.6)

xlabel('T_{peak}','FontName','Times New Roman','FontSize',fs,'Interpreter','tex')
ylabel('T_{night}','FontName','Times New Roman','FontSize',fs,'Interpreter','tex')
zlabel('UTCI_{peak}','FontName','Times New Roman','FontSize',fs,'Interpreter','tex')

ax = gca;
ax.FontName = 'Times New Roman';
ax.FontSize = fs;
ax.Box = 'on';
ax.LineWidth = 0.75;
grid on
ax.GridAlpha = 0.15;
view(135,25)

colormap(parula)
cb = colorbar;
cb.Label.String = 'GP predictive \sigma';
cb.Label.FontName = 'Times New Roman';
cb.Label.FontSize = fs;
cb.FontName = 'Times New Roman';

legend({'GP-sampled front (color = \sigma)', ...
    '144 samples',sprintf('True non-dominated (still, n=%d)', sum(trueStillNonDominated)), ...
        sprintf('True, now dominated (n=%d)', sum(~trueStillNonDominated))}, ...
    'FontName','Times New Roman','FontSize',fs-1,'Location','northoutside','Box','off')


%% ===================== HELPER FUNCTIONS =====================
function J = chebyshevCost(x, Mdl_night, Mdl_day, Mdl_utci, y_min, y_span, lambda, rho)
    y = [predict(Mdl_night,x), predict(Mdl_day,x), predict(Mdl_utci,x)];
    y_norm = (y - y_min) ./ y_span;
    weighted_diff = abs(y_norm) .* lambda;
    J = max(weighted_diff) + rho*sum(weighted_diff);
end

