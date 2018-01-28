function P = to_polytope(obj)
%TO_POLYTOPE Convert the given hyperrectangle to a polytope object.
%
%   P = to_polytope(H)
%       Convert hyperrectangle H to polytope P (array to array)
%
% (C) 2011 by Truong X. Nghiem (nghiem@seas.upenn.edu)

P = polytope;

if obj.dims == 0
    return;
end

H = [eye(obj.dims); -eye(obj.dims)];

n = size(obj.L, 2);
allP = cell(1, n);

K = [obj.H; -obj.L];

for kk = 1:n
    allP{kk} = polytope(H, K(:, kk), 1, 1);
end

P = horzcat(allP{:});

end
