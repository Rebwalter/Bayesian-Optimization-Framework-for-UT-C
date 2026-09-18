function mask = paretoFront(Y)
%PARETOFRONT Identify non-dominated (Pareto-optimal) rows of an objective matrix.
%   MASK = PARETOFRONT(Y) returns a logical column vector the same
%   height as Y, where MASK(i) is true if row Y(i,:) is not dominated
%   by any other row. All objectives in Y are assumed to be
%   minimized. Row i is dominated by row j if Y(j,:) <= Y(i,:) for
%   every column and Y(j,:) < Y(i,:) for at least one column.
%
%   Inputs:
%       Y - [n x p] matrix of p objective values for n design points
%
%   Outputs:
%       mask - [n x 1] logical vector, true for non-dominated rows

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
