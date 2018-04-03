function slipenchuk_problem5f
%SLIPENCHUK_PROBLEM5B
% **DISCLAIMER: This code is a modified version of 
% temple_abm_population_predator_prey.m located at 
% https://math.temple.edu/~seibold/teaching/2018_2100/
%
% ---Prompt---
% (b) Now change the speed of warriors of one army to 0.01, i.e., make the warriors 
% move slower.  Run the simulation, explain the difference in outcome, and 
% explain why this happens.
% (c) Change the probability of death to p death = min {0.01m,1}, and run 
% your code with both warrior speeds equal to 0.02 (case (a)), as well as 
% with warrior speeds from case (b).  Describe and explain the simulation
% results.
% (d) Change the probability of death to p death = min{0.1?m,1}, and run your 
% code with both warrior speeds equal to 0.02 (case (a)), as well as with warrior 
% speeds from case (b).  Describe and explain the simulation results.
% (e) Argue which real situations the three different laws of death 
% (in (a) vs. (c) vs. (d)) could model in terms of agent interactions.
% ---Implementation---
% Outcomes of simulation
% (b) Left Warrior speeds = .01, Right Warrior speeds = .02. Initially, I see that the 
% battle line moves to the left usually by a magnitude of 1. Since the left army is
% moving slower but their warrior spawning is unchanged, they can keep
% their territory from being invaded because their area is more dense. When
% the left army does go into the right side of the domain(past 5 x-axis)
% the right army is fast to take down their warriors, at first. This is 
% because the left army's slower motion prevents them from gaining the density they 
% need to maintain this new territory. Once their forces 'catch up' they
% can actually maintain the area. In this way the left army continues to
% gain territory in the right side of the domain.
%
% (c) Warrior speeds = .02. I noticed that warriors are able to spend
% relatively large amounts of time behind enemy lines. I believe this could
% simulate cavalry members charging past melee warriors and other enemy
% cavalry. The armies territories are fairly equal and each has their
% enemies relatively deep into their own territories. 
% (c) Left Warrior speeds = .01, Right Warrior speeds = .02. I noticed
% similiar results in the first part of (c) in that there are warriors, I
% think cavalry, that are able to rush through enemy lines on both armies.
% The battle line does not really change like in the first part. The left
% army isn't able to take territory despite its slow speed. I believe the
% changed probability of death is reason why. Despite their increased
% density of the left army, the death probability is not high enough for
% the left army to clear away the right army. 
%
% (d) Warrior speeds = .02. Here it looks like the armies are mostly
% composed of cavalry that travel deep into eachothers territories and
% fight their mini battles theire as opposed to a concentrated section like
% in part (b). The starting domains of each armies have larger densities of
% their own soldiers but the ratio becomes more even towards the center of
% the domain. This seems to be an exaggerated form of part (c) where only a
% few warriors, cavalry, could break through enemy lines. 
% (c) Left Warrior speeds = .01, Right Warrior speeds = .02. The right army
% is able to penetrate into left territory faster and deeper. Eventually
% the left army advances into the right armies domain but they cannot
% maintain numbers there. The right army has no problem maintaining forces
% in the left army's domain. 
%
% (e) I think situation (a) models regualr melee combat. I believe that
% situation (c) involves cavalry that go behind enemy lines and survive
% relatively longer. I believe that (d) could either represent simply more
% cavalry or it could represent a zoomed in perspective on the battle line
% where the armys meet. I imagine 2 barbarian forces running into each
% other and fighting in their seperate battles as opposed to a central one
% along a 'vertical' line. In general, it 
% seems that when the speed is set slower this could signify defensive
% tactics as opposed to when it is daster than could be offensive tactics. 
% 03/2018 by Matthew Slipenchuk
%            tuf91673@temple.edu

% Parameters
N = [100 100]; % initial number of warrior L and warrior R
ns = 10000; % number of steps
ax = [0 10 0 10]; % problem domain
v = [.01,.02]; % speeds of warrior L and warrior R (distance per step)
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

D{1} = ones(N(1), 1)*0; % initial angles of direction of warrior_l
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
    cell_of_warrior_l = floor(X{1})*[1;ncx]+1; % index of cell of warrior_l agent
    cell_of_warrior_r = floor(X{2})*[1;ncx]+1; % cell index warrior_r agent
    number_of_warrior_l = zeros(ncx*ncy,1); % initialize count per cell
    number_of_warrior_rs = number_of_warrior_l; % initialize count per cell
    
    % Number of warrior_l in each cell
    [cell_indices,~,agents_in_cell] = unique(cell_of_warrior_l);
    number_of_agents_in_cell = accumarray(agents_in_cell,1);
    number_of_warrior_l(cell_indices) = number_of_agents_in_cell;
    % Number of warrior_rs in each cell
    [cell_indices,~,agents_in_cell] = unique(cell_of_warrior_r);
    number_of_agents_in_cell = accumarray(agents_in_cell,1);
    number_of_warrior_rs(cell_indices) = number_of_agents_in_cell;

    % Death probabilities of each agent
    % part b
    %p_warrior_l_death = min((.001*(number_of_warrior_rs(cell_of_warrior_l).^2)),1);
    %p_warrior_r_death = min((.001*(number_of_warrior_l(cell_of_warrior_r).^2)),1);
    % part c
    %p_warrior_l_death = min((.001*(number_of_warrior_rs(cell_of_warrior_l))),1);
    %p_warrior_r_death = min((.001*(number_of_warrior_l(cell_of_warrior_r))),1);
    % part d
    p_warrior_l_death = min((.001*(number_of_warrior_rs(cell_of_warrior_l).^(1/2))),1);
    p_warrior_r_death = min((.001*(number_of_warrior_l(cell_of_warrior_r).^(1/2))),1);
    
    % Population dynamics of warrior L
    ind_reproduce = 1;
    X_new = [1*(rand(1)/2), ax(4)*(rand(1))]; % copy the reproducing agents
    D_new = rand(sum(ind_reproduce),1)*0; % random directions
    ind_survive = rand(size(D{1}))>=p_warrior_l_death; % survining agents
    X{1} = [X{1}(ind_survive,:);X_new]; % new warrior_l positions
    D{1} = [D{1}(ind_survive);D_new]; % new warrior_l directions

    % Population dynamics of warrior R
    ind_reproduce = 1;
    X_new = [9 + (10)*(rand(1)),ax(4)*(rand(1))];% copy the reproducing agents
    D_new = rand(sum(ind_reproduce),1)*pi; % random directions
    ind_survive = rand(size(D{2}))>=p_warrior_r_death; % survining agents
    X{2} = [X{2}(ind_survive,:);X_new]; % new warrior_r positions
    D{2} = [D{2}(ind_survive);D_new]; % new warrior_r directions
    
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
    title(sprintf('warrior l (%d agents) and warrior r (%d agents)',...
        size(X{1},1),size(X{2},1)))
    drawnow
end
