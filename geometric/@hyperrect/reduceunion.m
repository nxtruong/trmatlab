function obj = reduceunion(obj)
%REDUCEUNION Remove redundant hyper-rectangles in a union (array).
%      B = reduceUnion(A)
%
%   All redundant hyper-rectangles in A (i.e. ones that are contained in
%   others) are removed. B is the result.
%
% (C) 2011 by Truong X. Nghiem (nghiem@seas.upenn.edu)

% Remove thin sets
obj = removethin(obj);

ndims = obj.dims;
if ndims == 0, return; end

% Sort the sets by their volumes
[theVolumes, idx] = sort(prod(obj.H - obj.L), 'descend');
obj.L = obj.L(:, idx);
obj.H = obj.H(:, idx);

k = 1;
while k < size(obj.L, 2)
    % Iterates from the largest sets; obj.Array may change because sets may
    % be removed from it.
    
    % Algorithm: use bsxfun() to find all sets after k that are contained
    % in set k (by comparing l and h); then use find() to find the indices
    % of those; then added k because they are mapped to sets after k in the
    % original Array.
    % This code is about 10 times faster than using loops.
    idx = find(...
        all(bsxfun(@ge, obj.L(:, k+1:end), obj.L(:, k)) & ...
            bsxfun(@le, obj.H(:, k+1:end), obj.H(:, k))))...
        + k;
    obj.L(:, idx) = [];
    obj.H(:, idx) = [];
    theVolumes(idx) = [];
    k = k + 1;
end


% More rigorous reduction: we loop through sets from small to large, and
% remove sets that is covered by the other sets.
N = length(theVolumes);
if N < 2, return; end

toKeep = true(1, N); % The sets that are being kept (have not been removed)
for k = N:-1:1
    toCheck = toKeep;
    toCheck(k) = false;  % The sets that are to be checked to cover set k
    idx = find(toCheck);
    
    % If nothing to check then quit
    if isempty(idx), break; end
    
    % If the total volume of the other sets is smaller than the current
    % set, they cannot cover it.
    if sum(theVolumes(toCheck)) < theVolumes(k)
        continue;
    end    
    
    CL = obj.L(:, k);
    CH = obj.H(:, k);
    
    % Do not check sets that do not intersect the current set
    idx(any(bsxfun(@ge, obj.L(:, idx), CH) | bsxfun(@le, obj.H(:, idx), CL))) = [];
    

    % Check for each sets in the union
    for kk = idx
        [CL, CH] = regiondiff(CL, CH, obj.L(:, kk), obj.H(:, kk));
        
        % Remove thin sets from (CL, CH)
        if ~isempty(CL)
            thinidx = any(CL >= CH);
            CL(:, thinidx) = [];
            CH(:, thinidx) = [];
            
            % Check if (CL, CH) is empty then we can stop
            if isempty(CL)
                break;
            end
        else
            break;
        end
    end
    
    if isempty(CL)
        % If the set difference is empty, i.e. other sets cover set k, then
        % we remove set k
        toKeep = toCheck;
    end
end

% Remove those sets that are not kept
toKeep = ~toKeep;
obj.L(:, toKeep) = [];
obj.H(:, toKeep) = [];

end
