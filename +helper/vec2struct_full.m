function s = vec2struct_full(x, gpParamNames)
%VEC2STRUCT_FULL Convert a numeric design vector into the UT&C parameter struct.
%   S = VEC2STRUCT_FULL(X, GPPARAMNAMES) maps a 1 x d numeric row
%   vector X (in the order given by GPPARAMNAMES) back into a struct
%   with one field per parameter, as expected by UTC.m. Radius_tree is
%   not itself a GP input; it is a derived quantity, reconstructed
%   here from the crown-radius fraction (alpha) and canyon width.
%
%   Inputs:
%       x            - [1 x d] numeric design point
%       gpParamNames - {1 x d} cell array of parameter names, in the
%                      same order as x
%
%   Outputs:
%       s - struct with one scalar field per parameter in
%           gpParamNames, plus a derived Radius_tree field

for j = 1:numel(gpParamNames)
    s.(gpParamNames{j}) = x(j);
end
s.Radius_tree = s.alpha * s.Width_canyon;
end
