function slipenchuk_problem5a
%SLIPENCHUK_PROBLEM5A
% **DISCLAIMER: This code is a modified version of 
% temple_abm_population_predator_prey.m located at 
% https://math.temple.edu/~seibold/teaching/2018_2100/
%
% ---Prompt---
% (a) Modify the Matlab file temple_abm_population_predator_prey.m
% from the course website
% http://math.temple.edu/~seibold/teaching/2018_2100/
% into a simulation of an iron age battlefield, as follows:
%  (i)  Agents represent warriors.  The two species represent two armies.  
% One army starts in (x,y)?[0,1]×[0,10].  The other starts in (x,y)?[9,10]×[0,10].  
% Place 100 agents on each side initially.
%  (ii)  There is no reproduction of agents.  Instead, in each step, one new 
% warrior is spawned at a random position in the domains given in (i).
%  (iii)  The probability of death of a warrior during a simulation step is
% p death = min{0.001m2,1}, where m is the number of warriors of the opposing 
% army in the same square.
%  (iv)  Warriors  starting  on  the  left  have  a  bias  to  move  to  the  right, 
% and  vice  versa.   This  is  modeled  as follows.  In each step, the direction angle
% d is a warrior is changed to 0.99d+ 0.01d0, where d0= 0 for warriors starting 
% on the left, and d0=? for warriors starting on the right.
%  (v)  For  now,  let  the  magnitude  of  angle  change  for  each  army  be  0.1,  
% and  the  speeds  of  warriors  from either army be 0.02.  That means that 
% the battle is symmetric, both in terms of geometric setup and in terms of 
% properties of the warriors/armies. Run the code multiple times, and 
% describe and explain the outcome of the simulation.  Submit (i.e., email
% to the course instructor) your program under the filename
% yourfamilyname_problem5a.m
% ---Implementation---
% Outcomes of simulation
% I notice that the armies meet in the middle of the domain and one will
% usually begin gaining territory before falling back and vice versa. This
% motion can be called a wave like. I also noticed that sometimes some
% warriors will break through the enemies lines and go deep into their
% territory before dieing. This is due to the randomness added to the
% probability of death. 
%
% 03/2018 by Matthew Slipenchuk
%            tuf91673@temple.edu

% Parameters
N = [100 100]; % initial number of warrior L and warrior R
ns = 10000; % number of steps
ax = [0 10 0 10]; % problem domain
v = [.02,.02]; % speeds of warrior L and warrior R (distance per step)
a = [.1;.1]; % magnitude of angle change per step for warrior L and warrior R

% Initialization
ncx = (ax(2)-ax(1)); % number of cells in horizontal direction
ncy = (ax(4)-ax(3)); % number of cells in vertical direction
cx = ax(1):ax(2); 
cy = ax(3):ax(4); % cell boundaries


X{1} = [1*(rand(N(1),1)/2),... % initial positions
    ax(4)*(rand(N(1),1))]; % Warrior L (left side)

X{2} = [9 + (10)*(rand(N(2),1)),... % initial positions
    ax(4)*(rand(N(2),1))]; % of Warrior R (right side)

D{1} = ones(N(1), 1)*0; % initial angles of direction of prey
D{2} = ones(N(1), 1)*pi; % initial angles of direction of warrior R

for j = 1:ns % loop over steps
    % Update positions and directions
    for k = 1:2 % warrior L and warrior R
        X{k}(:,1) = X{k}(:,1)+v(k)*cos(.99*D{k});
        X{k}(:,2) = X{k}(:,2)+v(k)*sin(.99*D{k} + (.01*pi)); % position
        D{k} = D{k}+a(k)*randn(size(D{k})); % change direction of motion
        
        % Let agents bounce of walls
        ind = (X{k}(:,1)<ax(1)&cos(D{k})<0)|... % who is hitting a wall
            (X{k}(:,1)>ax(2)&cos(D{k})>0); % horizontally
        D{k}(ind) = pi-D{k}(ind); % reverse x-direction
        ind = (X{k}(:,2)<ax(3)&sin(D{k})<0)|... % who is hitting a wall
            (X{k}(:,2)>ax(4)&sin(D{k})>0); % vertically
        D{k}(ind) = -D{k}(ind); % reverse y-direction
        
        % Move agents that are outside of domain into the domain
        X{k}(:,1) = min(max(X{k}(:,1),ax(1)+1e-9),ax(2)-1e-9);
        X{k}(:,2) = min(max(X{k}(:,2),ax(3)+1e-9),ax(4)-1e-9);
    end
    
    % Determine population balance in each cell
    cell_of_warrior_l = floor(X{1})*[1;ncx]+1; % index of cell of prey agent
    cell_of_warrior_r = floor(X{2})*[1;ncx]+1; % cell index predator agent
    number_of_prey = zeros(ncx*ncy,1); % initialize count per cell
    number_of_predators = number_of_prey; % initialize count per cell
    
    % Number of prey in each cell
    [cell_indices,~,agents_in_cell] = unique(cell_of_warrior_l);
    number_of_agents_in_cell = accumarray(agents_in_cell,1);
    number_of_prey(cell_indices) = number_of_agents_in_cell;
    % Number of predators in each cell
    [cell_indices,~,agents_in_cell] = unique(cell_of_warrior_r);
    number_of_agents_in_cell = accumarray(agents_in_cell,1);
    number_of_predators(cell_indices) = number_of_agents_in_cell;

    % Death probabilities of each agent
    p_prey_death = min((.001*(number_of_predators(cell_of_warrior_l).^2)),1);
    p_predator_death = min((.001*(number_of_prey(cell_of_warrior_r).^2)),1);
    
    % Population dynamics of warrior L
    ind_reproduce = 1;
    X_new = [1*(rand(1)/2), ax(4)*(rand(1))]; % copy the reproducing agents
    D_new = rand(sum(ind_reproduce),1)*0; % random directions
    ind_survive = rand(size(D{1}))>=p_prey_death; % survining agents
    X{1} = [X{1}(ind_survive,:);X_new]; % new prey positions
    D{1} = [D{1}(ind_survive);D_new]; % new prey directions

    % Population dynamics of warrior R
    ind_reproduce = 1;
    X_new = [9 + (10)*(rand(1)),ax(4)*(rand(1))];% copy the reproducing agents
    D_new = rand(sum(ind_reproduce),1)*pi; % random directions
    ind_survive = rand(size(D{2}))>=p_predator_death; % survining agents
    X{2} = [X{2}(ind_survive,:);9 + (10)*(rand(1,1)),... % initial positions
    ax(4)*(rand(1,1))]; % new predator positions
    D{2} = [D{2}(ind_survive);D_new]; % new predator directions
    
    % Plotting
    clf
    plot([1;1]*cx,ax([3,4])'*(cx*0+1),'k-',... % draw boundaries
        ax([1,2])'*(cy*0+1),[1;1]*cy,'k-') % of cells
    hold on
    plot(X{1}(:,1),X{1}(:,2),'.','markersize',12,'color',[0 .8 0])
    plot(X{2}(:,1),X{2}(:,2),'.','markersize',12,'color',[1 0 0])
    hold off
    axis equal xy, axis(ax)
    xlabel('x'), ylabel('y') % axis labels
    title(sprintf('Prey (%d agents) and predators (%d agents)',...
        size(X{1},1),size(X{2},1)))
    drawnow
end
