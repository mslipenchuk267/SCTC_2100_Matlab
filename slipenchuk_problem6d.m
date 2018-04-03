function slipenchuk_problem6d
%SLIPENCHUK_PROBLEM6D
% **DISCLAIMER: This code is a modified version of 
% temple_abm_traffic_car_following.m located at 
% https://math.temple.edu/~seibold/teaching/2018_2100/
%
% ---Prompt---
% (d) Change your code from part (c) to have dt = 1e-2 and np = 10
% again, but keep the Runge-Kutte 4 time stepping.  With this code, 
% conduct a parameter study how the characteristics of the traffic wave depend
% on the strength of the optimal velocity term.  Specifically, vary the coefficient 
% (which is 0.5 in the original file) in front of the term from 0 to 2, 
% in steps of 0.05.  For each choice, obtain the maximum and minimum vehicle 
% velocity, umax and umin,  at  the  final  time.   Using  these  data,  plot
% umax and umin as  functions  of the parameter.  
% Describe your observations and explain the result.  
% Submit your code under the filename yourfamilyname_problem6d.m
% ---Implemenation---
% (d) The time it takes for 'choke points' to appear is much faster the
% smaller the coeffiecent is. Meaning fastest at 0 and slowest at 2. I am
% confused as to how to plot the min and max velocities as functions of the
% coefficent. 
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
dt = 1e-2; % time step (for integration)
np = 10; % number of time steps between two plotting events
V = @(d) 10*(tanh(d/2-2)+tanh(2))/(1+tanh(2)); % optimal velocity function
f = @(x) [x(n+1:2*n);... % ODE right hand side
    20*(x([n+2:2*n,n+1])-x(n+1:2*n))./([x(2:n);x(1)+L]-x(1:n)-lv).^2+...
    2*(V([x(2:n);x(1)+L]-x(1:n)-lv)-x(n+1:2*n))];

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
min_final_vel = min(x(n+1:2*n)) % Final minimum velocity
max_final_vel = max(x(n+1:2*n)) % Final maxiumum velocity