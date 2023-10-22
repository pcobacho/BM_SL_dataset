% *************************************************************************
% AUTHORS: Pablo Cobacho
% *************************************************************************
% GROUP: Lab 1.3.5., Communications and Signal Processing Lab (ComSP), 
% TELMA, ETSIT, University of Malaga
% *************************************************************************
% DESCRIPTION:
% This script simulates the P2 beam management procedure in a mresobile 
% communications scenario with hexagonal cells. In addition, it calculates
% the downlink SINR for the users contained in the cell under study 
% (refCellID).
% *************************************************************************

clear; clc; close all;
rng(211);   % Set RNG state for repeatability

% Add search path with the project files
% *************************************************************************
[~, oldPath] = addPaths();

% Scenario parameters
scen_center = [0,0];  % center of scenario
num_users = 500;      % number of UEs in reference cell
interSiteDist = 200;  % inter site distance
refCellID = 13;        % reference cell ID (cell under study)
numScat = 0;         % number of scatters in scenario

% Get base stations (gNB) coordinates
gNBpos = get_gNB_positions(scen_center,interSiteDist);

% Get the cencell centers of all cells in scenario
cellCenters = getCellCenter(gNBpos,interSiteDist);

% Get user positions in reference cell
[userPos,~] = get_hex_users_positions(num_users,...
    cellCenters(:,mod(refCellID-1,3)+1,ceil(refCellID/3)), interSiteDist);

% Get scatters position
scatPos = get_scatters_positions(numScat,scen_center,interSiteDist);

% Plot scenario
fillCell = 1;
show_scenario(gNBpos,interSiteDist,cellCenters,fillCell,userPos,scatPos)

txBeamID = zeros(1,num_users);
rxBeamID = zeros(1,num_users);
optRSRP = zeros(1,num_users);
for u=1:num_users
    disp(['UserID: ' num2str(u)])
    [optRSRP(1,u),txBeamID(1,u),rxBeamID(1,u)] = get_serving_gNB_power([gNBpos(:,5);0], ...
        refCellID,[userPos(:,u);0],[scatPos;zeros(1,length(scatPos))]);
end

IDs = sort(unique(txBeamID));
colors = lines(length(IDs)+1);
figure; hold on, grid on, axis equal
for i=IDs
    uIdx = find(txBeamID==i); %user Index
    plot(userPos(1,uIdx),userPos(2,uIdx),'Color',colors(i,:),'Marker','o','LineStyle', 'none')
end


% Restore search paths
% *************************************************************************
path(oldPath);