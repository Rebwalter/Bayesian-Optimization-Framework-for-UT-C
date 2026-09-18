function area = area2D_dominated(points, r2, r3)
%AREA2D_DOMINATED Dominated area in 2D; the inner integral of HYPERVOLUME3D.
%   AREA = AREA2D_DOMINATED(POINTS, R2, R3) computes the area dominated
%   by a 2D point set, bounded above by reference values R2 and R3.
%
%   Inputs:
%       points - [n x 2] point set (objectives 2 and 3 of a 3D slice)
%       r2, r3 - scalar reference bounds for each dimension
%
%   Outputs:
%       area - scalar dominated area

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
