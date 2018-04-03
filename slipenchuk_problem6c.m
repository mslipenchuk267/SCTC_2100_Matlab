function slipenchuk_problem6c
%SLIPENCHUK_PROBLEM6C
% **DISCLAIMER: This code is a modified version of 
% temple_abm_traffic_car_following.m located at 
% https://math.temple.edu/~seibold/teaching/2018_2100/
%
% ---Prompt---
% (c) For the previous code, change the time stepping from Euler’s method 
% to Runge-Kutta 4.  Submit your code under the filename 
% yourfamilyname_problem6c.m, and explain why the new time stepping fixes the
% problem encountered in (b).
% ---Implementation---
% (c) I believe that implementing rk4 removes the previous issue because it
% is able to manage the larger time steps because it takes 4 the derivative
% 4 times once at an inital and endpoint, and twice at the midpoint. I
% think this gives a more accurate result whereas euler's method cannot
% predict as accurately what the velocity and position should be.
%
%
% 03/2018 by Matthew Slipenchuk
%            tuf91673@temple.edu

% Parameters
n = 22; % number of vehicles
L = 230; % length of road
w = 2; % width of road (for plotting only)
lv = 4.5; % length of each vehicle
tf = 60*5; % final time
dt = 1e-1; % time step (for integration)
np = 1; % number of time steps between two plotting events
V = @(d) 10*(tanh(d/2-2)+tanh(2))/(1+tanh(2)); % optimal velocity function
f = @(x) [x(n+1:2*n);... % ODE right hand side
    20*(x([n+2:2*n,n+1])-x(n+1:2*n))./([x(2:n);x(1)+L]-x(1:n)-lv).^2+...
    .5*(V([x(2:n);x(1)+L]-x(1:n)-lv)-x(n+1:2*n))];

% Initialization
q = linspace(0,L,n+1)'; q = q(1:end-1); % initial positions of vehicles
v = (1:n)'; % initial velocities of vehicles
x = [q;v]; % initial state vector
nt = ceil(tf/dt/np)*np; % number of time steps (muliple of np)
dt = tf/nt; % actual time step
phi = linspace(0,2*pi,101); cx = cos(phi); cy = sin(phi); % circle
r = L/(2*pi); % radius of road

% Computation
t = 0; % initial time
for j = 0:nt/np % time loop
    % Computation
    for i = 1:(np*(j>0)) % Do np compute steps (first time: do nothing)
        % Don't update x each time, i'm just taking a bunch of euler
        % forward steps
        s1 = f(x); % first slope
        s2 = f(x+dt/2*s1); % second slope
        s3 = f(x+dt/2*s2); % third slope
        s4 = f(x+dt*s3); % fourth slope
        x = x+dt*(s1+2*s2+2*s3+s4)/6;
        t = t+dt; % advance time
    end
    % Plotting
    q = x(1:n); % vehicle positions
    qv = [q-lv,q]'*2*pi/L; % vehicle front and back position angle
    clf
    subplot(1,2,1) % plot vehicles on circular road
    plot((r+w)*cx,(r+w)*cy,'k-',(r-w)*cx,(r-w)*cy,'k-') % draw road
    hold on
    plot(r*cos(qv),r*sin(qv),'b-','linewidth',3)
    hold off
    axis equal tight
    title(sprintf('Car-following model at t=%0.0fs: positions',t))
    subplot(1,2,2) % plot velocities over positions on line
    plot(mod(q,L),x(n+1:2*n),'b.')
    axis([0 L 0 20])
    title(sprintf('Car-following model at t=%0.0fs: velocities',t))
    xlabel('vehicle front position [m]')
    ylabel('vehicle velocity [m/s]')
    drawnow
end
