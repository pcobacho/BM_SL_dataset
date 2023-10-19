% *************************************************************************
% AUTHORS: Pablo Cobacho
% *************************************************************************
% GROUP: Lab 1.3.5., Communications and Signal Processing Lab (ComSP), 
% TELMA, ETSIT, University of Malaga
% *************************************************************************
% DESCRIPTION: From the positions of the base stations, the centres of 
% each cell and each user in the scenario, this scrip displays a graph of 
% the scenario with all its elements.
% *************************************************************************

function show_scenario(gNBpos,interSiteDist,cellCenters,fillCell,userPos)

r = interSiteDist/3;
num_gNB = length(gNBpos);
numCells = num_gNB*3;

% Convert hexadecimal color code to RGB (range 0 to 1)
colors = [sscanf('FF8080', '%2x').'; sscanf('D2E0FB', '%2x').';...
    sscanf('CDFAD5', '%2x').';sscanf('FFCF96', '%2x').';...
    sscanf('F9F3CC', '%2x').';sscanf('BEADFA', '%2x').';...
    sscanf('EDB7ED', '%2x').'] / 255;

figure();
hold on;

% Plot gNBs
plot(gNBpos(1,:), gNBpos(2,:), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'w'); 

% Plot cells
angles = linspace(0, 2*pi, 7);
for i=1:numCells
    vertPos = [r*cos(angles)+cellCenters(1,i); r*sin(angles)+cellCenters(2,i)];
    if fillCell
        fill(vertPos(1,:), vertPos(2,:), colors(ceil(i/3),:));
    else
        plot(vertPos(1,:), vertPos(2,:),'k-')
    end
    text(cellCenters(1,i), cellCenters(2,i), num2str(i), 'FontSize', 12);
end

% Plot users
if nargin==5
    scatter(userPos(1,:),userPos(2,:),'r.')
end

axis equal;
title('Hexagonal Cells Scenario');
xlabel('X (m)');
ylabel('Y (m)');
grid on;

hold off;