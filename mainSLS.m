% *************************************************************************
% AUTHORS: Pablo Cobacho
% *************************************************************************
% GROUP: Lab 1.3.5., Communications and Signal Processing Lab (ComSP), 
% TELMA, ETSIT, University of Malaga
% *************************************************************************
% DESCRIPTION:
% This script simulates the P2 beam management procedure in a mobile 
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
scen_center = [0,0];
num_users = 200;
interSiteDist = 200;
refCellID = 15;

gNBpos = get_gNB_positions(scen_center,interSiteDist);

cellCenters = getCellCenter(gNBpos,interSiteDist);

[userPos,~] = get_hex_users_positions(num_users,...
    cellCenters(:,mod(refCellID-1,3)+1,ceil(refCellID/3)), interSiteDist);

fillCell = 0;
show_scenario(gNBpos,interSiteDist,cellCenters,fillCell,userPos)

% Restore search paths
% *************************************************************************
path(oldPath);