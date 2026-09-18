function [DOE,X,Params] = makeDOE_opt(m,rdmseed)

Params = helper.defineParameters_opt();
paramNames = fieldnames(Params);
nP = numel(paramNames);

rng(rdmseed)
X = helper.lhsdesignNEW(m, nP);

for i = 1:m
    for j = 1:nP
        p = paramNames{j};
        lb = Params.(p).lb;
        ub = Params.(p).ub;
        DOE(i).(p) = lb + X(i,j) * (ub - lb);
    end

    % --------------------------------------------------
    % Tree geometry correction
    % --------------------------------------------------

    Hc = DOE(i).Height_canyon;
    Wc = DOE(i).Width_canyon;
    Dt = 0.1;
    Ht = 7.5;
    % Maximum radius allowed from geometry
    Rt_max1 = Hc / 2;
    Rt_max2 = (Wc - 2*Dt) / 4;
    Rt_max3 = Ht;
    Rt_max_allowed = min([Rt_max1, Rt_max2,Rt_max3]);
    % Convert alpha to radius
    alpha = DOE(i).alpha;

    Rt = alpha * Wc;

    % Rescale if too large
    if Rt > Rt_max_allowed
        Rt = 0.95 * Rt_max_allowed;  % small safety margin
    end

    DOE(i).Radius_tree = Rt;

    % ---- Height constraint ----
   

    Ht_max_allowed = Hc - Rt;

    if Ht > Ht_max_allowed
        Ht = 0.95 * Ht_max_allowed;
    end

    if Ht <= Rt
        Ht = Rt + 0.1;
    end



end

end

