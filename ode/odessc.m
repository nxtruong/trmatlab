function [varargout] = odessc(solver, ssc, method, odefun, tspan, y0, varargin)
%ODESSC ODE solver with custom step-size control
%   This function is a wrapper for Matlab's ODE solver functions (ode*)
%   with custom step-size control. At each step, a custom step-size control
%   function (ssc) is called to determine the maximum allowable value for
%   the next step size. The exact value of the next step size is still
%   determined by the ODE solver, with the upper-bound given by ssc.
%
%   Syntax:
%       varargout = odessc_step(solver, ssc, odefun, tspan, y0,...)
%   
%   Inputs:
%       solver - the chosen ODE solver (ode*), a string or function handle
%       ssc - step-size control
%       method - method to ensure the step-size control
%       tspan - must be of the form [t0 tfinal] where t0 < tfinal
%       odefun, y0,... - same arguments as the standard ODE solver
%
%   Outputs: same as the standard ODE solver
%
%   Step-size Control: ssc can be either
%       + a function name or a handle to a function (ssc_step) which
%           receives 2 input arguments: t and y, and returns the maximum
%           allowable value for the next step size. Note that if scc_step
%           returns a value which is too small, the solver may be stuck. So
%           it's a good practice to cap the value with a lower bound (e.g.
%           0.001).
%       + a cell array of two function handles [ssc_step ssc_check]:
%           ssc_step is the same as described above.  ssc_check receives
%           3 inputs: t, y, and tau; it returns a logical value
%           (true/false, 1/0) which tells if the next step size tau is
%           acceptable.  ssc_check is used to avoid calling ssc_step
%           many times (which is usually slow), so ssc_check should be
%           significantly faster than ssc_step.
%
%   Method is a number:
%       0: the solver is called once to obtain the entire solution, then
%           more samples might be interpolated if necessary to ensure
%           step-size control; this is the fastest method but also the less
%           accurate.
%       1: the solver is called multiple times to obtain the solution; this
%           method may be much slower than the previous, but more accurate;
%           however, it may generate too many samples.
%       2: similar to method "1" but generate fewer samples by removing
%           unnecessary ones; it is slower than method "1" with the same
%           accuracy.
%
% (C) 2009 by Truong Nghiem (nghiem@seas.upenn.edu)

% Check and process the input arguments
if nargin < 5, error('Not enough input arguments.'); end

if ischar(solver)
    solver = str2func(solver);
elseif ~isa(solver, 'function_handle')
    error('solver must be either a function name or a function handle.');
end

has_ssccheck = false;
if iscell(ssc)
    ssc_step = ssc{1};
    ssc_check = ssc{2};
    has_ssccheck = true;
elseif isa(ssc, 'function_handle')
    ssc_step = ssc;
elseif ischar(ssc)
    ssc_step = str2func(ssc);
else
    error('ssc must be either a function name or a function handle or a cell array of two function handles.');
end

if ~isvector(tspan) || length(tspan) ~= 2 || tspan(1) >= tspan(2)
    error('tspan must be of the form [t0 tfinal] where t0 < tfinal.');
end


% Solve the ODE
if method == 0
    x = []; y = [];  % Results
    
    % Call the ODE solver for the entire horizon
    sol = solver(odefun, tspan, y0, varargin{:});
    
    % Interpolate the samples to ensure the step-size control
    k = 1;  % the current position in sol.x
    nx = length(sol.x);
    
    while k < nx
        % Interpolate between sol.x(k) and sol.x(k+1)
        x = [x sol.x(k)];
        nextx = sol.x(k+1);
        cury = sol.y(:,k);
        y = [y cury];
        if ~has_ssccheck || ~ssc_check(x(end), cury, nextx-x(end))
            % Next step size is not OK --> split it
            nextstep = x(end) + ssc_step(x(end), cury);  % next step by ssc
            
            % While not exceed the next simulated step --> interpolate samples
            while nextstep < sol.x(k+1)
                % Add it to the result
                x = [x nextstep];
                cury = deval(sol, nextstep);
                y = [y cury];
                
                if has_ssccheck && ssc_check(nextstep, cury, nextx-nextstep)
                    % Next step size is OK --> quit
                    break;
                end
                nextstep = nextstep + ssc_step(nextstep, cury); % next step by ssc
            end
        end
        
        k = k + 1;
    end
    
    % Add final point
    x = [x sol.x(end)];
    y = [y sol.y(:, end)];
    
    % Events
    if isfield(sol, 'xe')
        te = sol.xe;
        ye = sol.ye;
        ie = sol.ie;
    else
        [te, ye, ie] = deal([]);
    end
    
elseif method == 1 || method == 2
    method1 = (method == 1);
    
    cury = y0;
    tcurrent = tspan(1);
    tfinal = tspan(2);
    
    x = tcurrent; y = y0;  % Results
    te = []; ye = []; ie = [];

    runningsolver = true;
    
    while runningsolver
        nextstep = tcurrent + ssc_step(tcurrent, cury);  % Next step by ssc
        
        if nextstep >= tfinal
            runningsolver = false;
            nextstep = tfinal;
        end  % cannot exceed tfinal
        
        % Call the solver
        sol = solver(odefun, [tcurrent nextstep], cury, varargin{:});
        
        % Extract results
        % The first point is tcurrent, so we remove it
        
        % if the solver stopped because of a terminal event, stop it
        if sol.x(end) < nextstep, runningsolver = false; end
        
        if method1
            x = [x sol.x(2:end)];
            y = [y sol.y(:, 2:end)];
        else
            if length(sol.x) == 2 || ~runningsolver || ...
                    (~isempty(sol.xe) && sol.xe(end) == nextstep)
                % Only 1 new point
                % Or gonna stop now
                % Or if an event happened at nextstep --> save all
                x = [x sol.x(2:end)];
                y = [y sol.y(:, 2:end)];
            else
                % length > 2 --> save until the point just before the last
                x = [x sol.x(2:end-1)];
                y = [y sol.y(:, 2:end-1)];
            end
        end
        tcurrent = x(end);
        cury = y(:, end);
        
        % Save the events
        te = [te sol.xe];
        ye = [ye sol.ye];
        ie = [ie sol.ie];
    end
else
    error('Unrecognized method.');
end


% Return results to varargout
if nargout == 1
    % Return sol structure
    [sol.x, sol.y, sol.xe, sol.ye, sol.ie] = deal(x, y, te, ye, ie);
    varargout{1} = sol;
elseif nargout > 1
    outputs = {x.', y.', te.', ye.', ie.'};
    [varargout{1:nargout}] = outputs{1:nargout}; % populate the outputs
end

end