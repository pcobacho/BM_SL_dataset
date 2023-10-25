% *************************************************************************
% AUTHORS: Pablo Cobacho
% *************************************************************************
% GROUP: Lab 1.3.5., Communications and Signal Processing Lab (ComSP), 
% TELMA, ETSIT, University of Malaga
% *************************************************************************
% DESCRIPTION: this script calculates the position of a given number of 
% users (num_users) located in the cell under study.
% *************************************************************************

function [userPos,vertPos] = get_hex_users_positions(num_users,cellCenter,interSiteDist,hUE)

r = interSiteDist/3; % hexagonal cell radius
angles = linspace(0, 2*pi, 7); % angles of the vertices of the hexagonal cell
vertPos = [cellCenter(1)+r*cos(angles); cellCenter(2)+r*sin(angles)]; % vertice coordinates

% generate random points
userPos(1,:) = cellCenter(1)+(2*rand(1,num_users)-1) * r;
userPos(2,:) = cellCenter(2)+(2*rand(1,num_users)-1) * r;
userPos(3,:) = hUE*ones(1,num_users);

insidePos = inpolygon(userPos(1,:), userPos(2,:), vertPos(1,:), vertPos(2,:));
numOutsidePos = sum(~insidePos);
while numOutsidePos > 0
    userPos(1,~insidePos) = cellCenter(1) + (2*rand(1,numOutsidePos)-1) * r;
    userPos(2,~insidePos) = cellCenter(2) + (2*rand(1,numOutsidePos)-1) * r;
    userPos(3,~insidePos) = hUE;
    insidePos = inpolygon(userPos(1,:), userPos(2,:), vertPos(1,:), vertPos(2,:));
    numOutsidePos = sum(~insidePos);
end