%% INITIAL LHS DESIGN AND UT&C MODEL EVALUATION

% This script generates the initial Latin Hypercube Sampling (LHS)
% design, evaluates the UT&C model for all parameter combinations,
% and calculates the three objective functions used for
% multi-objective optimisation:
%
%   1. Nighttime minimum air temperature
%   2. Daytime peak air temperature
%   3. Peak UTCI
%
% The resulting dataset provides the initial training data for the
% subsequent Bayesian optimisation.
%
% Required files:
%   - MeteoData_DWD_hourly_RadPart.mat
%   - helper.makeDOE_opt.m
%   - UTC.m
%
% Outputs:
%   - DOE
%   - Params
%   - y_night
%   - y_peak
%   - y_UTCI_peak
%   - MeteoData
%   - Results2m
%   - UTCI

%% Configuration ==========================================================

meteoFile = '+data/+input/MeteoData_DWD_hourly_RadPart.mat';

meteoRows = 1009:1105; % Heat wave 2021

% Study site
Name_Site = 'HST_BEM';

% Initial design
nInitialSamples = 50;
randomSeed = 12;

% Time window used to determine the nighttime minimum temperature
% following the hottest day.
nightWindow = 40:60;

%% Load meteorological data ===============================================

load(meteoFile, 'MeteoData_DWD');

MeteoData = MeteoData_DWD(meteoRows, :);

%% Generate initial LHS design ============================================

[DOE, ~, Params] = helper.makeDOE_opt( ...
    nInitialSamples, randomSeed);

%% Start parallel pool ====================================================

if isempty(gcp('nocreate'))
    parpool('local');
end

poolObj = gcp();
nWorkers = poolObj.NumWorkers;

fprintf('Using %d parallel workers\n', nWorkers);

%% Split DOE into chunks ==================================================

chunkEdges = round(linspace(0, nInitialSamples, nWorkers+1));

DOE_chunks = cell(nWorkers,1);
n_chunks   = zeros(nWorkers,1);

for c = 1:nWorkers
    idxRange = (chunkEdges(c)+1):chunkEdges(c+1);

    DOE_chunks{c} = DOE(idxRange);
    n_chunks(c)   = numel(idxRange);
end

%% Parallel UT&C evaluation ===============================================

T2m_chunks  = cell(nWorkers,1);
UTCI_chunks = cell(nWorkers,1);

parfor c = 1:nWorkers

    DOE_c = DOE_chunks{c};
    n_c   = n_chunks(c);

    % Run UT&C for this chunk
    [Results2m_c, UTCI_c] = UTC( ...
        DOE_c, n_c, MeteoData, Name_Site);

    % Convert air temperature from Kelvin to degrees Celsius
    Results2m_c.T2m = Results2m_c.T2m - 273.15;

    % Store only the required outputs
    T2m_chunks{c}  = Results2m_c.T2m;
    UTCI_chunks{c} = UTCI_c;

end

%% Reassemble results in original DOE order ===============================

Results2m.T2m = cat(3, T2m_chunks{:});
UTCI          = cat(3, UTCI_chunks{:});

%% Sanity checks ==========================================================

assert(size(Results2m.T2m,3) == nInitialSamples, ...
    'Reassembled T2m has %d slices, expected %d.', ...
    size(Results2m.T2m,3), nInitialSamples);

assert(size(UTCI,3) == nInitialSamples, ...
    'Reassembled UTCI has %d slices, expected %d.', ...
    size(UTCI,3), nInitialSamples);

%% Calculate objective functions ==========================================

% -------------------------------------------------------------------------
% Objective 1: Peak daytime air temperature
% -------------------------------------------------------------------------

T_peak = max(Results2m.T2m, [], 1);
y_peak = reshape(T_peak, [], 1);

% -------------------------------------------------------------------------
% Objective 2: Minimum nighttime air temperature
% -------------------------------------------------------------------------

T_window = Results2m.T2m(nightWindow, :, :);
T_night = min(T_window, [], 1);
y_night = reshape(T_night, [], 1);

% -------------------------------------------------------------------------
% Objective 3: Peak UTCI
% -------------------------------------------------------------------------

UTCI_peak = max(UTCI, [], 1);
y_UTCI_peak = reshape(UTCI_peak, [], 1);

%% Combine objective values ===============================================

% Objective order:
%   [T_night, T_peak, UTCI_peak]

y_all = [y_night, y_peak, y_UTCI_peak];

%% Save initial dataset ===================================================

save('+data/+results/Initial_LHS_Dataset.mat', ...
    'DOE', ...
    'Params', ...
    'MeteoData', ...
    'Name_Site', ...
    'nInitialSamples', ...
    'randomSeed', ...
    'y_night', ...
    'y_peak', ...
    'y_UTCI_peak', ...
    'y_all', ...
    'Results2m', ...
    'UTCI');