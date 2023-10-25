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
num_users = 300;      % number of UEs in reference cell
hUE = 0;              % UEs height
num_gNB = 7;          % number of base stations (gNBs)
h_gNB = 25;            % gNB height
num_cells = num_gNB*3; % number of cells (sectors)
interSiteDist = 200;  % inter site distance
refCellID = 1;        % reference cell ID (cell under study)
numScat = 0;         % number of scatters in scenario

% Get base stations (gNB) coordinates
gNBpos = get_gNB_positions(scen_center,interSiteDist,h_gNB);

% Get the cencell centers of all cells in scenario
cellCenters = getCellCenter(gNBpos,interSiteDist);

% Get user positions in reference cell
[userPos,~] = get_hex_users_positions(num_users,...
    cellCenters(:,mod(refCellID-1,3)+1,ceil(refCellID/3)), interSiteDist,hUE);

% Get scatters position
scatPos = get_scatters_positions(numScat,scen_center,interSiteDist);

% Plot scenario
fillCell = 0;
show_scenario(gNBpos,interSiteDist,cellCenters,fillCell,userPos,scatPos)

% Obtains the optimal beam pair RSRP each user of the reference cell
txBeamID = zeros(1,num_users);
rxBeamID = zeros(1,num_users);
optRSRP = zeros(1,num_users);
fprintf('<strong>Obtaining optimal beam pair RSRP of each user:</strong>\n');
for u=1:num_users    
    disp(['UserID: ' num2str(u)])
    
    current_gNBpos = gNBpos(:,ceil(refCellID/3));
    [~,optRSRP(1,u),txBeamID(1,u),rxBeamID(1,u)] = get_serving_gNB_power(current_gNBpos, ...
        refCellID,userPos(:,u),scatPos);
end

IDs = sort(unique(txBeamID));
colors = lines(length(IDs)+1);
cIdx=1; % color index
figure; hold on, grid on, axis equal
for i=IDs
    uIdx = find(txBeamID==i); %user Index
    plot(userPos(1,uIdx),userPos(2,uIdx),'Color',colors(cIdx,:),'Marker','o','LineStyle', 'none')
    cIdx = cIdx+1;
end
% 
% % Calculates interfering RSRP
% num_beams = 8;
% beamID = randi(num_beams,1,num_cells);
% rsrp = zeros(num_users,num_cells);
% fprintf('<strong>Obtaining interfering RSRP of each user:</strong>\n');
% for u=1:num_users
%     disp(['UserID: ' num2str(u)])
%     for c=1:num_cells
%         current_gNBpos = gNBpos(:,ceil(c/3));
%         rsrp(u,c) = get_interfering_gNB_rsrp(current_gNBpos,c,userPos(:,u),...
%             scatPos,beamID(c));
%     end
% end


% Restore search paths
% *************************************************************************
path(oldPath);