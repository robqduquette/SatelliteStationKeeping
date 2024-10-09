clear
close all

%% Station keeping example problem

% constants
earth.r = 6378e3; % m
earth.mu = 0.39860e15; % m3/s2

target.rho = 10 * 1e-3 * 1e-9; % kg/m3
target.alt_des = 400e3; % m
target.alt_band = 1e1; % m, the tolerance for the desired altitude
target.r = target.alt_des + earth.r;
target.period = 2*pi* sqrt(target.r^3/earth.mu);
target.r_max = target.r + target.alt_band;
target.r_min = target.r - target.alt_band;

sc.thrust = 0.017; % N
sc.mass = 1000; % kg
sc.a_thrust = sc.thrust/sc.mass;
sc.CD = 0.3;
sc.cross_section = 10; %m2
sc.a_drag = @(r) 0.5 * target.rho * earth.mu * sc.CD * sc.cross_section / (sc.mass * r);

problem.sc = sc;
problem.target = target;
problem.body = earth;

%% simulation setup
% y = (r, along_track, ap, thrusting)
y0 = [target.r, 0, -1]; % Start with coasting state
periods = 10;
t_end = target.period * periods;

% Set 0up time and event detection
options = odeset('Events', @(t,y) target_event(t,y,problem));

% Initialize time and solution arrays
t_total = [];
y_total = [];
ye_total = [];
time_elapsed = 0;
% Loop to handle state updates and continue integration after events
while true
    % Integrate the system until an event occurs
    [t, y, te, ye, ie] = ode45(@(t,y) sat_dynamics(t,y,problem), [0 t_end], y0, options);
    
    t_real = t + time_elapsed;

    % Append results to total solution
    t_total = [t_total; t_real];
    y_total = [y_total; y];
    ye_total = [ye_total; ye];
    
    % If no event occurs, break the loop
    if isempty(te)
        break;
    end
    
    % adjust time to real time
    time_elapsed = time_elapsed + t(end);

    % Update state based on event type
    r_event = ye(1);
    
    if r_event <= problem.target.r_min
        disp(['Boosting at time ', num2str(te/3600), ' hours']);
        % Boosting state: set thrust to active (1) and reset altitude slightly above min
        y0 = [ye(1)+0.001, ye(2), 1];
    elseif r_event >= problem.target.r_max
        disp(['Coasting at time ', num2str(te/3600), ' hours']);
        % Coasting state: set thrust to inactive (-1) and reset altitude slightly below max
        y0 = [ye(1) - 0.001, ye(2), -1];
    end
    
    % Adjust time span to simulate the remaining time
    t_end = t_end - te;
end

% Plot the results
r = y_total(:,1);
along_track = y_total(:,2);
x = r - target.r;

plot(along_track, x, '-b')
hold on
scatter(ye_total(:,2), ye_total(:,1)-target.r,'or')
scatter(along_track(1), x(1), 'og')

grid on
xlabel('Along Track Error (m)')
ylabel('SMA Error (m)')
title('Relative Satellite motion during station keeping')

% time graph
figure
plot(t_total/3600, y_total(:,1) - target.r)
hold on
plot(t_total/3600, y_total(:,2))
xlabel('time (hours)')
ylabel('position error (meters)')
legend('radial error','along track error')
grid on
title(['Satellite Position error over ', num2str(periods),' orbits'])


%% functions
% y = (r, along_track, ap, thrusting)
function dydt = sat_dynamics(t, y, problem)
    thrust = y(3) > 0;
    r = y(1);    
    mu = problem.body.mu;
    ap = problem.sc.a_thrust * thrust - problem.sc.a_drag(r);

    drdt = 2 * ap * r^1.5 / sqrt(mu);
    dtrackdt = sqrt(mu / r) - sqrt(mu / problem.target.r);

    dydt = [drdt; dtrackdt; 0];
end

function [value, isterminal, direction] = target_event(t, y, problem)
    r = y(1);
    value = [r - problem.target.r_min; r - problem.target.r_max];
    isterminal = [1; 1]; % Stop the integration at event
    direction = [0; 0]; % Detect both directions (crossing either boundary)
end
