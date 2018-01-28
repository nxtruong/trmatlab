function obj = combine(obj)
%COMBINE Combine hyperrectangles in a union.
%      B = combine(A)
%
%   Combine hyperrectangles in an array A, if possible, to reduce the
%   number of them.  B is the result.
%
%   This function can be very heavy, so only use it when necessary.
%
% (C) 2011 by Truong X. Nghiem (nghiem@seas.upenn.edu)

% Remove thin sets
obj = removethin(obj);

ndims = obj.dims;
if ndims == 0, return; end

N = size(obj.L, 2);

combined = true;
while combined   % repeat until no combination was made
    if N < 2, return; end

    combined = false;
    
    k = 1;
    while k < N
        % Find sets that are the same as k except for at most 1 dimension
        idx = find(sum(bsxfun(@eq, obj.L(:, k), obj.L(:, k+1:end)) &...
                       bsxfun(@eq, obj.H(:, k), obj.H(:, k+1:end))...
                   ) >= ndims-1, 1);
        if ~isempty(idx)
            % Combine them
            combined = true;
            
            idx = idx + k;
            obj.L(:, k) = min(obj.L(:, k), obj.L(:, idx));
            obj.H(:, k) = max(obj.H(:, k), obj.H(:, idx));
            
            obj.L(:, idx) = [];
            obj.H(:, idx) = [];
            
            N = N - 1;
        end
        
        k = k + 1;
    end
end

end
