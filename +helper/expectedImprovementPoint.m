function EI = expectedImprovementPoint(x, GP, f_best, varNames)
%EXPECTEDIMPROVEMENTPOINT Expected Improvement (EI) at a single design point.
%   EI = EXPECTEDIMPROVEMENTPOINT(X, GP, F_BEST, VARNAMES) evaluates
%   the Expected Improvement acquisition function at design point X
%   (a 1 x d row vector) for a minimization problem, using the
%   predictive mean/variance of a trained fitrgp model.
%
%   Inputs:
%       x        - [1 x d] design point
%       GP       - trained fitrgp surrogate of the scalarized cost
%       f_best   - scalar, best (minimum) observed cost so far
%       varNames - variable names; kept for interface parity with
%                  callers that build x from a table, unused here
%                  since predict() takes a plain numeric row vector
%
%   Outputs:
%       EI - scalar expected improvement (0 if predictive sigma is
%            effectively zero)

[mu, sigma] = predict(GP, x);
sigma = max(sigma, 1e-6);

improvement = f_best - mu;
Z = improvement ./ sigma;
EI = improvement .* normcdf(Z) + sigma .* normpdf(Z);

if sigma < 1e-8
    EI = 0;
end
end
