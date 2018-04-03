function slipenchuk_problem4a
%SLIPENCHUK_PROBLEM4A
% **DISCLAIMER: This code is a modified version of 
% temple_abm_population_migrate_mate_and_age.m located at 
% https://math.temple.edu/~seibold/teaching/2018_2100/
%
% ---Prompt---
% (a) Modify the Matlab files temple_abm_population_migrate_mate_and_age.m
% from the course website http://math.temple.edu/~seibold/teaching/2018_2100/
% as follows: (a) increase the speed of the agents to 0.02;(b) place initially 
% 500 agents uniformly in the left half domain, (x, y)?[0,5]×[0,10];(c) change 
% the code so that offspring producing cells are ones that contain 7, 8, 9, or 10 agents; 
% (d) change the code so that (in each step) an offspring producing cell produces 
% a new agent with probability 1/4.  Run your code multiple times and explain 
% your observations.  Submit (i.e., email to the course instructor) your program
% under the filename yourfamilyname_problem4a.m
%
% (b) Change the probability of a cell producing a new agent to 1/5.  
% Describe and explain what the model produces, and why.
% (c) Leave the probability of a cell producing a new agent at 1/5, but now 
% let a cell be offspring producing if it contains 7, 8, 9, 10, 11, 12, 13, 14, 
% or 15 agents.  Describe and explain what the model produces, and why.
% ---Implementation---
% --Observations--
% The agents sturuggle to cross over to the right hand side of the domain
% but eventually works its way over.
% This stems from the fact that it now takes 7-10 agents in a cell to create
% a new agent, the starting location of the agents, and the low probability
% of reproduction.
%
% The increased of amount of agents required for reproduction affects the
% population growth's spread. I observed that the agents have a hard time 
% spreading to the right side of the x domain. There are natrually less 
% agents there, and this makes it harder for a cell to contain 7-10 agents.
% I also noticed that on the left hand side of the x domain, this
% requirement keeps the agents population from 'over mating'. This will
% sometimes lead to some holes or pockets appearing on the left hand side.
%
% The starting location of the agents effects this spread as well. If the
% agents were starting in the middle of the x and y axis, they would spread
% uniformly in each direction. However, the agents start in the left side
% of the grid. The cells closer to the left have less issue reproducing
% than those on the right. 
%
% The 25% probability of reproduction further decreases the rate that the
% population can spread right. The few cells on the right that do have 7-10
% cells now only have a 25% chance to reproduce. 
%
% So the increased amount of agents required to reproduce, the starting
% location, and the reproduction probability all slow, but do not ultimately stop,
% the spread of the agents to the right.
%
% 02/2016 by Matthew Slipenchuk
%            tuf91673@temple.edu

% Parameters
N = 500; % initial umber of agents
ns = 5000; % number of steps
ax = [0 10 0 10]; % problem domain
v = .02; % speed of agents (distance traveled per step)
max_age = 100; % each agents lives 100 steps
p = .20; % Probability for mating .25 for part a, .20 for part b

% Initialization
cx = ax(1):ax(2); cy = ax(3):ax(4); % cell boundaries
X = [ax(1)+(ax(2)-ax(1))/2*rand(N,1),... % initial positions
    ax(4)*rand(N,1)]; % of agents (left half of domain)
D = rand(N,1)*2*pi; % initial angles of direction of agent
A = randi(max_age,N,1); % age of agents

for j = 1:ns % loop over steps
    % Update positions and directions
    X(:,1) = X(:,1)+v*cos(D); X(:,2) = X(:,2)+v*sin(D); % move agents
    D = D+0.1*randn(size(D)); % change direction of motion
    A = A+1; % increase age of each agent by 1
    
    % Let agents bounce of walls
    ind = (X(:,1)<ax(1)&cos(D)<0)|... % who is hitting a wall
        (X(:,1)>ax(2)&cos(D)>0); % horizontally
    D(ind) = pi-D(ind); % reverse x-direction
    ind = (X(:,2)<ax(3)&sin(D)<0)|... % who is hitting a wall
        (X(:,2)>ax(4)&sin(D)>0); % vertically
    D(ind) = -D(ind); % reverse y-direction
    X(:,1) = min(max(X(:,1),ax(1)),ax(2)); % move agents outside of
    X(:,2) = min(max(X(:,2),ax(3)),ax(4)); % domain back onto boundary
    
    % Mating: if more than one agent in same cell, create new agent
    p = rand(1);
    if p <= .25 % Probability offspring ready cells reproduce
        [occupied_cells,~,agents_in_cell] = unique(floor(X),'rows');
        cells_with_multiple_agents = occupied_cells(...
            accumarray(agents_in_cell,1)>=7 & ...
            accumarray(agents_in_cell,1)<=10,:);
        N_new = size(cells_with_multiple_agents,1); % number of new agents
        X = [X;cells_with_multiple_agents+.5]; % new agents in cell centers
        D = [D;rand(N_new,1)*2*pi]; % random direction
        A = [A;zeros(N_new,1)]; % age of new agents is 0
    end
    % Remove agents who are too old
    ind = A<=max_age; % agents who do not die
    X = X(ind,:); D = D(ind); A = A(ind); % keep those
    
    % Plotting
    clf
    plot([1;1]*cx,ax([3,4])'*(cx*0+1),'k-',... % draw boundaries
        ax([1,2])'*(cy*0+1),[1;1]*cy,'k-') % of cells
    hold on
    ind = A<=max_age/3; % young agents
    plot(X(ind,1),X(ind,2),'.','markersize',12,'color',[0 .8 0])
    ind = max_age/3<A&A<=max_age*2/3; % middle aged agents
    plot(X(ind,1),X(ind,2),'.','markersize',12,'color',[.8 .8 0])
    ind = max_age*2/3<A; % old agents
    plot(X(ind,1),X(ind,2),'.','markersize',12,'color',[1 0 0])
    hold off
    axis equal xy, axis(ax)
    xlabel('x'), ylabel('y') % axis labels
    title(sprintf('Migrating, mating, and aging animals (%d agents)',...
        length(A)))
    drawnow
    if ~numel(A), break, end % stop of no agents left
end
