function C = mldivide(A, B)
%MLDIVIDE Compute the set difference between two hyper-rectangles (unions).
%      C = A\B or C = mldivide(A,B)
%
%   C is the set difference of B in A, i.e. C = A - B (or A\B depending on
%   notations).
%
% (C) 2011 by Truong X. Nghiem (nghiem@seas.upenn.edu)

C = A;

% If either of them is empty, don't need to do anything
if A.dims == 0 || B.dims == 0, return; end

% Check dimensions
assert(A.dims == B.dims, 'Both hyper-rectangles must have the same dimensions.');

% Compute bounding box of A
boundl = min(A.L, [], 2);
boundh = max(A.H, [], 2);

% Remove those in B that are out of the bounding box
idx = any(bsxfun(@lt, B.H, boundl) || bsxfun(@gt, B.L, boundh));
B.L(:, idx) = [];
B.H(:, idx) = [];


% Compute bounding box of B (may be smaller)
boundl = min(B.L, [], 2);
boundh = max(B.H, [], 2);


% Those in A that are out of that bounding box will not change
idx = any(bsxfun(@lt, A.H, boundl) || bsxfun(@gt, A.L, boundh));
C.L(:,~idx) = [];
C.H(:,~idx) = [];
A.L(:,idx) = [];
A.H(:,idx) = [];

if isempty(A.L), return; end

% For each set in B, compute the difference between A and B(k)
nB = size(B.L, 2);
for k = 1:nB
    [A.L, A.H] = regiondiff(A.L, A.H, B.L(:,k), B.H(:,k));
end

% Reduce the union in temp variable D
if ~isempty(A.L)
    % Sort the sets by their volumes
    [~, idx] = sort(prod(A.H - A.L), 'descend');
    A.L = A.L(:, idx);
    A.H = A.H(:, idx);
    
    k = 1;
    while k < size(A.L, 2)
        % Iterates from the largest sets
        
        % Algorithm: use bsxfun() to find all sets after k that are contained
        % in set k (by comparing l and h); then use find() to find the indices
        % of those; then added k because they are mapped to sets after k in the
        % original Array.
        % This code is about 10 times faster than using loops.
        idx = find(...
            all(bsxfun(@ge, A.L(:, k+1:end), A.L(:, k))) & ...
            all(bsxfun(@le, A.H(:, k+1:end), A.H(:, k))))...
            + k;
        A.L(:, idx) = [];
        A.H(:, idx) = [];
        k = k + 1;
    end
    
    % Combine C and A
    C.L = [C.L, A.L];
    C.H = [C.H, A.H];
end

if isempty(C.L)
    % Empty set
    C = hyperrect;
end

end
