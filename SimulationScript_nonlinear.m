function SimulationScript()
    %% Initialize the simulation variables
    t0 = 0; % initial time
    dt = 0.1; % time step
    tf = 10.0; % final time
        
    % Set the control input function
    u = @(t, x) sin(t);    
    
    % Set the starting point
    x0 = [0;0;0;0];
    
    %% Simulate and plot the system using ode
    % Simulate the system
    [tvec, xvec] = matlabOde45(x0, t0, dt, tf, u);
    uvec = getControlVector(tvec, xvec, u);
    
    % Plot the resulting states
    figure;
    plotResults(tvec, xvec, uvec, 'b');
    
    %% Simulate and plot the system using Euler (or other method)
    % Simulate the system
    [tvec, xvec] = eulerIntegration(x0, t0, dt, tf, u);
    uvec = getControlVector(tvec, xvec, u);
    
    % Plot the results
    plotResults(tvec, xvec, uvec, 'r:');
    
    
end

function u_vec = getControlVector(tvec, xvec, u)
%getControlVector calculate the control vector over the specified time
%interval and state
%
% Inputs:
%   tvec: 1xm vector of time inputs
%   xvec: nxm matrix of states
%   u: function handle that takes time and state as inputs and outputs
%   the control input

    len = size(tvec, 2);
    u_vec = zeros(1, len);
    for k = 1:len
        u_vec(:,k) = u(tvec(k), xvec(:,k));
    end

end

function [tvec, xvec] = matlabOde45(x0, t0, dt, tf, u)
    %MatlabOde45 uses ODE 45 to simulate the state starting at x0 from time
    % t0 to tf
    %
    % Inputs:
    %   x0: nx1 initial state
    %   t0: scalar - initial time
    %   dt: scalar - time increment
    %   tf: scalar - final time
    %   u: function handle that takes time and state as inputs and outputs
    %   the control input
    %
    % Outputs:
    %   tvec: 1xm vector of times associated with the states
    %   xvec: nxm matrix of states where each column is a state at the
    %   associated time in tvec
    
    % Initialize the time
    t = t0:dt:tf;
    
    % Simulate the output
    [tvec xvec] = ode45(@(t,x) f(t,x,u(t,x)), t, x0);
    
    % Transpose the outputs to get in the correct form
    tvec = tvec';
    xvec = xvec';    
end

function [tvec, xvec] = eulerIntegration(x0, t0, dt, tf, u)
    %eulerIntegration uses eulerIntegration to simulate the state starting at x0 from time
    % t0 to tf
    %
    % Inputs:
    %   x0: nx1 initial state
    %   t0: scalar - initial time
    %   dt: scalar - time increment
    %   tf: scalar - final time
    %   u: function handle that takes time and state as inputs and outputs
    %   the control input
    %
    % Outputs:
    %   tvec: 1xm vector of times associated with the states
    %   xvec: nxm matrix of states where each column is a state at the
    %   associated time in tvec
    
    % Initialize values
    tvec = t0:dt:tf;
    len = length(tvec);
    xvec = zeros(size(x0,1), len);
    xvec(:,1) = x0;
    
    % Simulate forward in time
    % Write the euler simulation code here    
end

function zdot = f(t, z, u)
    %f calculates the state dynamics using the current time, state, and
    %control input
    %
    % Inputs:
    %   t: current time
    %   x: current state
    %   u: current control input
    %
    % Ouputs:
    %   xdot: time derivative of x(t)
    z2 = z(2);
    M = 0.5; % Mass of cart
    m = 0.2; % mass of pendulum
    b = 0.1; % coefficient of friction for cart
    l = 0.3; % length to pendulum center of mass
    I = 0.006; % mass moment of inertia of the pendulum
    g = 9.8; % Gravity constant

    denom = (I*(M+m) + M*m*l^2);
    A22 = (-(I+m*l^2)*b)/denom;
    A23 = (m^2*g*l^2)/denom;
    A42 = (-m*l*b)/denom;
    A43 = (m*g*l*(M+m))/denom;
    u2 = (I + m*l^2)/denom;
    u4 = (m*l)/denom;

    % Linear system matrices
    A = [0 1   0   0; 
         0 A22 A23 0; 
         0 0   0   1;
         0 A42 A43 0];
    B = [0; u2; 0; u4];

    % LTI equation
    zdot = A*z + B*u;
end

function plotResults(tvec, xvec, uvec, color)

    % Plot variables
    fontsize = 18;
    linewidth = 2;
    
    % Plot the resulting states
    subplot(5,1,1); hold on;
    plot(tvec, xvec(1,:), color, 'linewidth', linewidth);
    ylabel('x_1(t)', 'fontsize', fontsize);
    
    subplot(5,1,2); hold on;
    plot(tvec, xvec(2,:), color, 'linewidth', linewidth);
    ylabel('x_2(t)', 'fontsize', fontsize);
    
    subplot(5,1,3); hold on;
    plot(tvec, xvec(3,:), color, 'linewidth', linewidth);
    ylabel('x_3(t)', 'fontsize', fontsize);
    
    subplot(5,1,4); hold on;
    plot(tvec, xvec(4,:), color, 'linewidth', linewidth);
    ylabel('x_4(t)', 'fontsize', fontsize);
  
    subplot(5,1,5); hold on;
    plot(tvec, uvec, color, 'linewidth', linewidth);
    ylabel('u(t)', 'fontsize', fontsize);
    xlabel('time (s)', 'fontsize', fontsize);

end

