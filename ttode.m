function [varargout] = ttode(solver, events, odefun, tspan, y0, varargin)
%TTODE Solve ODE with time-based events
% Solve a system of differential equations with time-based and state-based
% events.  Time-based events (time-triggered) are events that happen at
% pre-defined time instants, for example periodic events.  State-based
% events happen when the state of the system satisfies certain conditions.
% State-based events are detected by zero-crossings (or roots) of the state
% trajectory.  This function is useful for simulating trajectories of
% switched systems, in which the dynamics of the system changes after every
% event.
%
% Syntax:
%   [TOUT, YOUT, TE, YE, IE] = TTODE(SOLVER, EVENTS, ODEFUN, TSPAN, Y0...)
%
% SOLVER is a function handle to a built-in ODE solver, e.g. @ode45, which
% will be used to solve the ODE.
%
% EVENTS is a structure with the following fields:
%   .ttevents is a cell array that specifies the time-based events. Each
%       cell is either a vector of at least 2 numbers or a function handle.
%       A vector [T, ofs1, ofs2, ..., ofsk] of (k+1) elements specifies k
%       periodic time events: ofs1 + jT, ofs2 + jT,... where j = 0,1,2,..
%       and T > 0 is the time period.
%       A function handle fevent(t, y) returns the next time event te > t,
%       where y is the current state.
%   .sevents specifies state-based events and is exactly the same as the
%       event function used in ODESET.  If the solver terminates because
%       one of the terminating state-based events has happened, this
%       function also terminates.  However, if a time-based event also
%       happened at exactly the same time as the terminating state-based
%       event did, this function will not terminate (i.e. time-based events
%       have higher priority than state-based events).
%   .zenoEps and .zenoN are used to detect Zeno phenomenons so that the
%       solver does not run forever (literally).  If there are zenoN > 0
%       consecutive events that are less than zenoEps time units apart from
%       the next, then the solver stops and returns -1 in IE(end).  If
%       zenoN <= 0, zenoEps <= 0, or these fields are omitted, no Zeno
%       detection will be performed.
%
%   In the future, callback functions may be implemented for the events.
%   For now, the state trajectory is continuous, thus this function cannot
%   simulate a bouncing ball (it is not designed to simulate hybrid systems
%   with state resets).
%
% All other input arguments are exactly the same as the built-in ODE
% solvers in Matlab.  Note that the Events field in the OPTIONS structure
% (set by ODESET) may be overwritten if EVENTS.sevents is defined.  So you
% should define only one of them, either in EVENTS.sevents (recommended) or
% in OPTIONS.
%
% All output arguments are the same as the built-in ODE solvers, except
% that it does not support the output solution structure (to be used with
% DEVAL).
%
% See also: ODE45, ODESET, ODEGET
% 
% (C) 2011 by Truong X. Nghiem (truong DOT nghiem AT gmail)
%                              (OR nghiem AT seas DOT upenn DOT edu)

% Process and Check input arguments
error(nargchk(5, inf, nargin));
error(nargchk(0, 5, nargout));

assert(isa(solver, 'function_handle') && ~isempty(strmatch('ode', func2str(solver))),...
    'SOLVER must be a function handle to one of the built-in ODE solver.');

assert(isstruct(events), 'EVENTS must be a structure; refer to the help.');

assert(numel(tspan) == 2 && tspan(2) > tspan(1),...
    'TSPAN must be of the form [T0, TFINAL] where TFINAL > T0.');

T0 = tspan(1);
TFINAL = tspan(2);

hasTTEvents = isfield(events, 'ttevents');
if hasTTEvents
    assert(iscell(events.ttevents), 'EVENTS.TTEVENTS must be a cell array.');
    assert(all(cellfun(...
        @(c) isvector(c) || isa(c, 'function_handle'), events.ttevents)),...
        'A cell in EVENTS.TTEVENTS must be either a vector or a function handle.');
    nTTEvents = numel(events.ttevents);
end

hasZeno = all(isfield(events, {'zenoEps', 'zenoN'})) && ...
    events.zenoEps > 0 && events.zenoN > 0;

if isfield(events, 'sevents')
    if ~isempty(varargin)
        % The OPTIONS structure exists
        options = odeset(varargin{1}, 'Events', events.sevents);
    else
        options = odeset('Events', events.sevents);
    end
    varargin{1} = options;
end


% Initialize

T = T0;
Y = y0(:)';

if nargout > 0
    tempout = {T0, Y, [], [], []};
    [varargout{1:nargout}] = tempout{1:nargout};
end

if nargout < 2
    tempout = cell(1, 2);  % at least get the time instants & states
else
    tempout = cell(1, nargout);
end


if hasZeno
    nZeno = 0;   % Number of possible Zeno events
    tZeno = T0;  % Time instant of the previous event
end

% Solve the ODE
while T < TFINAL
    % Calculate the next time event
    nextEvent = TFINAL;
    if hasTTEvents
        for k = 1:nTTEvents
            if isvector(events.ttevents{k})
                period = events.ttevents{k}(1);
                offsets =  events.ttevents{k}(2:end);
                
                % Check two steps to eliminate numerical errors
                ksteps = max(0, floor((T - offsets)/period + 1));
                allTs = ksteps * period + offsets;
                allTs(allTs <= T) = [];  % Remove erroneous entries
                nextT = min([allTs, ksteps * period + period + offsets]);
            else
                nextT = max(T, events.ttevents{k}(T, Y));
            end
            if nextT < nextEvent
                nextEvent = nextT;
            end
        end
    end
    
    % Simulate the system from T to nextEvent
    [tempout{:}] = solver(odefun, [T, nextEvent], Y, varargin{:});
    
    % Extract states and times
    T = tempout{1}(end);
    Y = tempout{2}(end,:)';
    
    % Save the outputs
    if nargout > 0
        % Save TOUT and YOUT
        for k = 1:min(nargout,2)
            varargout{k} = [varargout{k}; tempout{k}(2:end, :)];  % Cut the first row
        end
    end
    
    if nargout > 2
        % Save event results
        for k = 3:nargout
            varargout{k} = [varargout{k}; tempout{k}];
        end
    end
    
    % If the end time is < nextEvent, i.e. some terminating state event
    % happened, then this function should terminate
    if T < nextEvent
        break;
    end
    
    % Check for Zeno
    if hasZeno
        if T - tZeno < events.zenoEps
            nZeno = nZeno + 1;
        else
            nZeno = 0;
        end
        if nZeno >= events.zenoN
            if nargout >= 5
                varargout{5}(end + 1) = -1;
            end
            warning('TTODE:Zeno', 'Zeno phenomenon happens at time t = %g.', T);
            break;
        end
    end
end

% DONE

end