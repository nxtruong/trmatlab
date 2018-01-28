function obj = merge(obj, Options)
%MERGE Merge hyperrects in a union to reduce their numbers.
%      B = merge(A, Options)
%
%   There are two methods:
%   1) By expansion: each hyperrect in the union will be expanded in
%   different ratios. The largest ratio with that the expanded hyperrect is
%   still inside the union will be used. Then all other hyperrects that lie
%   inside the newly expanded one will be removed.
%   2) By using polytope/merge: the hyperrects are converted into
%   polytopes, which are then merged by calling polytope/merge, then the
%   resulted polytopes are converted back to hyperrects.
%
%   Options is a struct:
%   Options.method = 1 if method 1, 2 if method 2. Default: 1.
%   Options.ratioeps = the accuracy of method 1; a number between 0 and 1.
%                       Default 0.01.
%   If method 2 is used, Options may contain other fields for
%   polytope/merge.  See that function for details.
%
% See also: HYPERRECT/REDUCEUNION, POLYTOPE/MERGE
%
% (C) 2012 by Truong X. Nghiem (nghiem@seas.upenn.edu)

error(nargchk(1, 2, nargin));

% Remove thin sets
obj = removethin(obj);

ndims = obj.dims;
if size(obj.H, 2) < 2 || ndims == 0, return; end

defaultEps = 0.01;

if nargin < 2
    Options = struct('method', 1, 'ratioeps', defaultEps);
else
    assert(isstruct(Options), 'Options must be a structure.');
    
    if ~isfield(Options, 'method')
        Options.method = 1;
    else
        assert(Options.method == 1 || Options.method == 2, ...
            'Invalid method.');
    end
    
    if Options.method == 1
        if ~isfield(Options, 'ratioeps')
            Options.ratioeps = defaultEps;
        else
            assert(isscalar(Options.ratioeps) && ...
                Options.ratioeps > 0 && Options.ratioeps < 1, ...
                'Invalid Options.ratioeps.');
        end
    end
end

if Options.method == 2
    P = merge(to_polytope(obj), Options);
    
    % Convert P back to obj
    [H, K] = double(P);
    if iscell(H)
        obj.H = zeros(ndims, length(H));
        obj.L = obj.H;
        
        for kk = 1:length(H)
            [I, ~] = find(H{kk} > 0);
            obj.H(:,kk) = K{kk}(I);

            [I, ~] = find(H{kk} < 0);
            obj.L(:,kk) = -K{kk}(I);
        end
    else
        [I, ~] = find(H > 0);
        obj.H = K(I);
        
        [I, ~] = find(H < 0);
        obj.L = -K(I);
    end
    
    return;
end


% Sort the sets by their volumes
[theVolumes, idx] = sort(prod(obj.H - obj.L), 'descend');
obj.L = obj.L(:, idx);
obj.H = obj.H(:, idx);

% Bounding box
BL = min(obj.L, [], 2);
BH = max(obj.H, [], 2);

% The trivial case: if obj is exactly a hyperrect (its bounding box)
if rawsubset(BL, BH, obj.L, obj.H, false, false, false)
    obj.L = BL;
    obj.H = BH;
    return;
end

k = 1;
while k < size(obj.L, 2)
    % Iterates from the largest sets; obj.Array may change because sets may
    % be removed from it.
    
    % max lengths to expand: for lower bound for upper bound.
    maxL = obj.L(:,k) - BL;
    maxH = BH - obj.H(:,k);
    
    % Find the ratio with that the expanded hyperrect is still inside obj
    foundRat = 0;
    maxRat = 1;

    while maxRat - foundRat > Options.ratioeps
        rat = (foundRat + maxRat)/2;
        
        XL = obj.L(:,k) - rat * maxL;
        XH = obj.H(:,k) + rat * maxH;
        
        % Check if X is a subset of obj
        % Don't need to sort X and obj
        if rawsubset(XL, XH, obj.L, obj.H, false, false, false)
            foundRat = rat;
        else
            maxRat = rat;
        end
    end
    
    % Only if this hyperrect can be expanded, we remove all other sets
    % contained in this set
    if foundRat
        % Algorithm: use bsxfun() to find all sets after k that are contained
        % in set k (by comparing l and h); then use find() to find the indices
        % of those; then added k because they are mapped to sets after k in the
        % original Array.
        % This code is about 10 times faster than using loops.
        XL = obj.L(:,k) - foundRat * maxL;
        XH = obj.H(:,k) + foundRat * maxH;
        
        idx = find(...
                   all(bsxfun(@ge, obj.L, XL) & ...
                       bsxfun(@le, obj.H, XH)));

        idx(idx == k) = [];  % Remove k from idx
        
        if ~isempty(idx)
            % Accept XL, XH
            obj.L(:,k) = XL;
            obj.H(:,k) = XH;
            
            % Adjust volume of k
            theVolumes(k) = prod(XH - XL);
            
            obj.L(:, idx) = [];
            obj.H(:, idx) = [];
            theVolumes(idx) = [];
            
            % Adjust k if any sets before k are removed
            k = k - sum(idx < k);
        end
    end
    
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
