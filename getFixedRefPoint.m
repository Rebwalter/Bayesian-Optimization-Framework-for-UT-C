function r = getFixedRefPoint()
%GETFIXEDREFPOINT Fixed 3-objective hypervolume reference point.
%   R = GETFIXEDREFPOINT() returns the reference point [T_night_worst,
%   T_peak_worst, UTCI_peak_worst] used for EVERY hypervolume
%   calculation in this project (LHS baseline, each of the 10
%   independent BO runs, and the dense GP-reconstructed fronts), so
%   that all reported hypervolume values sit on a single common scale
%   and are directly comparable across scripts and runs.
%
%   Defined once here, rather than re-typed as a literal in each
%   script, so a future change only has to be made in one place.
%
%   Outputs:
%       r - [1 x 3] reference point, in objective order
%           [T_night, T_peak, UTCI_peak]

r = [20.5279, 35.1913, 41.8699];
end
