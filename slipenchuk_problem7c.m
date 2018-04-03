function slipenchuk_problem7c
%SLIPENCHUK_PROBLEM7C
% **DISCLAIMER: This code is a modified version of 
% temple_abm_cellular_game_of_life_file_ic.m located at 
% https://math.temple.edu/~seibold/teaching/2018_2100/
%
% ---Prompt---
% Use the Matlab file temple_abm_cellular_game_of_life_file_ic.m
% from the course website http://math.temple.edu/~seibold/teaching/2018_2100/
% to produce three interesting programs, named yourfamilyname_problem7a.m,
% yourfamilyname_problem7b.m, and yourfamilyname_problem7c.m,
% to be submitted as usual (i.e., email to the course instructor).
% Each  of  the  three  files  should  run  an  animation  of  Conway�s  Game  of  Life,  with  at  least  two  interesting
% objects interacting with each other in interesting ways, in the spirit of 
% the provided example of a glider gun shooting at a blinker ship.
% Take the structures from an online resource like the website
% http://www.radicaleye.com/lifepage/glossary.html (or other resources).
% Make sure to submit the corresponding *.txt files of the objects with your 
% codes.
% ---Implementation---
%   Implementation of Conway's Game of Life - Mayhem in the Kitchen. 
%   "A cell remains alive if it has 2 or 3 live neighbors, and a dead cell 
%   becomes alive if it has exactly 3 neighbors. Otherwise cells die or
%   remain dead." 
%   A group of teenagers set off a fuse igniting a pressure cooker located 
%   in a friends kitchen. 
%   I used leightweight_fuse.txt from
%   http://www.radicaleye.com/lifepage/glossary.html.
%   I used pressure_cooker.txt from
%   http://www.conwaylife.com/ref/lexicon/lex_p.htm#puffertrain.
%   These files should be included in the working directory:
%   lightweight_fuse.txt
%   pressure_cooker.txt
%
% 04/2018 by Matthew Slipenchuk
%            tuf91673@temple.edu

% Parameters
n = [250,100]; % number of cells per dimension
filename{1} = 'lightspeed_wire.txt'; % Fuse 
pos{1} = [110,45]; % position of object loaded from file (top corner)
filename{2} = 'pressure_cooker.txt'; % Pressure Cooker
pos{2} = [97,45]; % position of object loaded from file (top corner)
% Construct initial configuration
F = zeros(n); % initialize empty array
for j = 1:length(pos) % loop over objects
    S = load_data_file(filename{j}); % road data file
    F(pos{j}(1)+(1:size(S,1)),... % add object
        pos{j}(2)+(1:size(S,2))) = S; % to 2d array
end

% Initialization
n = size(F); % new size of domain (if necessary)
shl1 = [n(1),1:n(1)-1]; shr1 = [2:n(1),1]; % shift index vectors in dim 1
shl2 = [n(2),1:n(2)-1]; shr2 = [2:n(2),1]; % shift index vectors in dim 2

% Computation
clf
for j = 0:3000 % time loop
    % Plotting
    imagesc(~F') % plot array of cells
    axis equal tight, caxis([0,1]), colormap(gray)
    title(sprintf('Mild Act of Arson after %d steps',j))
    pause(1e-2+(j==0)) % wait a bit (and a bit more initially)
    % Update rule
    neighbors = F(shl1,shl2)+F(shl1,:)+F(shl1,shr2)+... % number
        F(:,shl2)+F(:,shr2)+... % of neighbors
        F(shr1,shl2)+F(shr1,:)+F(shr1,shr2); % for each cell
    F = (F&(neighbors==2|neighbors==3))|... % cell remains alive
        (~F&neighbors==3); % or cell becomes alive
end

%========================================================================

function S = load_data_file(filename)
% Read the data file filename and assign the information in 2d array S.
fid = fopen(filename,'r'); % open data file for reading
data = double(fscanf(fid,'%c',inf)); % read data as single string
fclose(fid); % close file
data = [10,data,10]; % add space in beginning and end of string
ind_data = data==42|data==46; % indices where actual cell data is stored
i_start = find(diff(ind_data)==1); % beginning index of row
i_end = find(diff(ind_data)==-1); % beginning index of row
row_length = i_end(1)-i_start(1); % length of each row in data
data = data(ind_data); % remove non-cell information
data = data==42; % make data logical for cell states
S = reshape(data,row_length,[]); % object stored as 2d array
