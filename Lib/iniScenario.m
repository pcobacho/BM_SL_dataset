function [gNBpos,cellCenters,userPos,scatPos] = iniScenario(prm)

% Get base stations (gNB) coordinates
gNBpos = get_gNB_positions(prm.scen_center,prm.interSiteDist,prm.h_gNB);

% Get the cencell centers of all cells in scenario
cellCenters = getCellCenter(gNBpos,prm.interSiteDist);

% Get user positions in reference cell
[userPos,~] = get_hex_users_positions(prm.num_users,...
    cellCenters(:,mod(prm.refCellID-1,3)+1,ceil(prm.refCellID/3)), ...
    prm.interSiteDist,prm.hUE);

% Get scatters position
scatPos = get_scatters_positions(prm.numScat,prm.scen_center, ...
    prm.interSiteDist);

if prm.showFigures
    % Plot scenario
    show_scenario(gNBpos,prm.interSiteDist,cellCenters,prm.fillCell,userPos,scatPos)
end