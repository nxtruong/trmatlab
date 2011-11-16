function display(obj)
%DISPLAY  Display a hyperrect object.
%
% (C) 2011 by Truong X. Nghiem (nghiem@seas.upenn.edu)

if obj.dims == 0
    disp('An empty hyper-rectangle.');
else
    N = length(obj);
    if N > 1
        fprintf('A union of %d hyper-rectangles in %d-D.\n', N, obj.dims);
    elseif obj.dims < 20
        % If dimension is not too large, display the full data
        fprintf('A hyper-rectangle in %d-D with intervals:\n', obj.dims);
        disp([obj.L, obj.H]);
    else
        fprintf('A hyper-rectangle in %d-D.\n', obj.dims);
    end
end

end
