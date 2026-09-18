function HV = hypervolume3D(P, r)
%HYPERVOLUME3D Compute the dominated hypervolume of a 3-objective Pareto set.
%   HV = HYPERVOLUME3D(P, R) computes the hypervolume of the region
%   dominated by the non-dominated point set P, bounded above by the
%   fixed reference point R. Uses the HSO (Hypervolume by Slicing
%   Objectives) algorithm: slice along objective 1, then integrate the
%   dominated area in the remaining two objectives at each slice.
%   All objectives are assumed to be minimized (lower is better).
%
%   Inputs:
%       P - [n x 3] non-dominated objective points
%       r - [1 x 3] reference point, dominated by every point in P
%
%   Outputs:
%       HV - scalar dominated hypervolume
%
%   Note: R should be held FIXED across all runs/comparisons that will
%   be compared against each other, so hypervolume values are on a
%   common scale (see ref_point_fixed in run_LHS_baseline_and_GP.m).

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
    if width <= 0
        continue
    end
    active = P(1:k, 2:3);
    area = helper.area2D_dominated(active, r(2), r(3));
    HV = HV + width * area;
end
end
