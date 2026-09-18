function J = chebyshevCost(x, Mdl_night, Mdl_day, Mdl_utci, y_min, y_span, lambda, rho)
%CHEBYSHEVCOST Augmented Tchebycheff scalarization, evaluated on GP surrogate predictions.
%   J = CHEBYSHEVCOST(X, MDL_NIGHT, MDL_DAY, MDL_UTCI, Y_MIN, Y_SPAN,
%   LAMBDA, RHO) evaluates the augmented Tchebycheff scalarization (as
%   used in ParEGO) at design point X, using the three trained GP
%   surrogate models to predict [T_night, T_peak, UTCI_peak]. This is
%   the objective minimized by fmincon when locating a single point on
%   the approximate Pareto front for a given weight vector LAMBDA.
%
%   Inputs:
%       x                   - [1 x d] design point
%       Mdl_night/day/utci  - trained fitrgp surrogate models
%       y_min, y_span       - [1 x 3] normalization constants (min, range),
%                             computed from the LHS baseline objectives
%       lambda              - [1 x 3] Sobol-sampled simplex weight vector
%       rho                 - scalar augmentation coefficient
%
%   Outputs:
%       J - scalar scalarized cost

y = [predict(Mdl_night,x), predict(Mdl_day,x), predict(Mdl_utci,x)];
y_norm = (y - y_min) ./ y_span;
weighted_diff = abs(y_norm) .* lambda;
J = max(weighted_diff) + rho*sum(weighted_diff);
end
