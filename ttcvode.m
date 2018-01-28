function [varargout] = ttcvode(events, odefun, tspan, y0, options, LMM, NLS, last)
%TTCVODE Solve ODE with time-based events using CVODE.
% Solve a system of differential equations with time-based events.
% Time-based events (time-triggered) are events that happen at pre-defined
% time instants, for example periodic events.  This function is useful for
% simulating trajectories of switched systems, in which the dynamics of the
% system changes after every event.
%
% Syntax:
% [TOUT, YOUT] = TTODE(EVENTS, ODEFUN, TSPAN, Y0, OPTIONS, LMM, NLS, LAST)
%
% This function uses CVODES, thus it requires sundialsTB - the Matlab
% interface of the Sundials library - to be installed.  See 
% http://computation.llnl.gov/casc/sundials for more information.
%
% EVENTS is a cell array that specifies the time-based events. Each cell is
%   either a vector of at least 2 numbers or a function handle.  A vector
%   [T, ofs1, ofs2, ..., ofsk] of (k+1) elements specifies k periodic time
%   events: ofs1 + jT, ofs2 + jT,... where j = 0,1,2,.. and T > 0 is the
%   time period.  A function handle fevent(t, y) returns the next time
%   event te > t, where y is the current state.
%
% ODEFUN is a function [dy, flag] = ODEFUN(t, y) which returns the
%   right-hand side of the ODE: y' = ODEFUN(t, y) and a flag (0 if OK).  It
%   may also have a third input argument which is the user data.  See
%   CVRhsFn for more information.
%
% TSPAN is of the form [t0 tfinal] where t0 < tfinal specifies the time
%   span of the simulation.
%
% Y0 is the initial condition (vector).
%
% OPTIONS is the option structure, which can be set by CVodeSetOptions.  It
%   can be used to specify the Linear Solver, Jacobian function, etc.  It
%   is optional.
%
% LMM and NLS are the solvers used by CVODES.  See CVodeInit for more
% information.  The defaults are LMM = 'BDF' and NLS = 'Newton'.
%
% LAST is true/1 then only the last solution is returned in TOUT and YOUT.
%   Otherwise, all solutions are returned.  Default: false.
%
% TOUT is a vector of solution time instants.
% YOUT is a matrix whose columns are the solutions at each time instant in
%   TOUT.  Note that this is different from the convention used by Matlab's
%   ODE solvers (ode45,...) and ttode, where the solutions are in the rows
%   of YOUT.
%
%
% See also: CVODEINIT, CVODESETOPTIONS, CVRHSFN
% 
% (C) 2011 by Truong X. Nghiem (truong DOT nghiem AT gmail)
%                              (OR nghiem AT seas DOT upenn DOT edu)

% HISTORY
%   2012-07-25  Truong fixed bug with ttevents containing function handles.

% Process and Check input arguments
error(nargchk(4, inf, nargin));
error(nargchk(0, 2, nargout));

assert(iscell(events), 'EVENTS must be a cell array; refer to the help.');
assert(all(cellfun(...
        @(c) (isnumeric(c) && isvector(c)) || isa(c, 'function_handle'), events)),...
        'A cell in EVENTS must be either a vector or a function handle.');
nEvents = numel(events);
    
assert(numel(tspan) == 2 && tspan(2) > tspan(1),...
    'TSPAN must be of the form [T0, TFINAL] where TFINAL > T0.');

T0 = tspan(1);
TFINAL = tspan(2);

if ~exist('LMM', 'var')
    LMM = 'BDF';
end

if ~exist('NLS', 'var')
    NLS = 'Newton';
end

if ~exist('options', 'var') || isempty(options)
   options = CVodeSetOptions('LMM', LMM, 'NonlinearSolver', NLS); 
end

if ~exist('last', 'var')
    last = false;
end

% Initialize

T = T0;
Y = y0(:);

% Reserve memory in blocks to improve performance
% For large systems, this helps improve the performance a lot.
if ~last
    MemBlk = min(floor(1024*1024/8/length(y0)), 256);
    NPoints = 1;  % number of points in the result, used to free excess memory
    curMem = MemBlk;
    
    if nargout >= 1
        varargout{1} = repmat(T0, MemBlk, 1);
    end
    
    if nargout >= 2
        varargout{2} = repmat(Y, 1, MemBlk);
    end
end

% Initialize CVODES
CVodeInit(odefun, LMM, NLS, T0, Y, options);


% Solve the ODE
status = 1;  % The status of the solver, 1 means 'compute next event'
while T < TFINAL
    if status == 1
        % Calculate the next time event
        nextEvent = TFINAL;
        
        for k = 1:nEvents
            if isnumeric(events{k})
                period = events{k}(1);
                offsets =  events{k}(2:end);
                
                % Check two steps to eliminate numerical errors
                ksteps = max(0, floor((T - offsets)/period + 1));
                allTs = ksteps * period + offsets;
                allTs(allTs <= T + eps) = [];  % Remove erroneous entries
                nextT = min([allTs, ksteps * period + period + offsets]);
            else
                nextT = max(T, events{k}(T, Y));
            end
            if nextT < nextEvent
                nextEvent = nextT;
            end
        end
        
        % Set new stop time
        CVodeSet('StopTime', nextEvent);
    end
    
    % Simulate the system for one step
    [status, T, Y] = CVode(TFINAL, 'OneStep');
    
    % If save all solutions
    if ~last
        % Increase memory if necessary
        if NPoints == curMem
            if nargout >= 1
                varargout{1}(end+1:end+MemBlk) = 0;
            end
            
            if nargout >= 2
                varargout{2}(:, end+1:end+MemBlk) = 0;
            end
            
            curMem = curMem + MemBlk;
        end
        
        % Save the outputs
        NPoints = NPoints + 1;
        
        if nargout >= 1
            varargout{1}(NPoints) = T;
            if nargout >= 2
                varargout{2}(:, NPoints) = Y;
            end
        end
    end
end

CVodeFree;

if last
    if nargout >= 1
        varargout{1} = T;
        if nargout >= 2
            varargout{2} = Y;
        end
    end    
else
% Free excess memory
if NPoints < curMem && nargout >= 1
    varargout{1}(NPoints+1:end) = [];
    if nargout >= 2
        varargout{2}(:, NPoints+1:end) = [];
    end
end
end

% DONE

end