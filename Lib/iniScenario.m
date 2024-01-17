function [gNBpos,cellCenters,userPos,scatPos] = iniScenario(prm)

% Get base stations (gNB) coordinates
gNBpos = get_gNB_positions(prm.scen_center,prm.interSiteDist,prm.h_gNB);

% Get the cencell centers of all cells in scenario
cellCenters = getCellCenter(gNBpos,prm.interSiteDist);

% Get user positions in reference cell
[userPos,~] = get_hex_users_positions(prm.num_users,...
    cellCenters(:,mod(prm.refCellID-1,3)+1,ceil(prm.refCellID/3)), ...
    prm.interSiteDist,prm.hUE);

r = 50;
theta =  linspace(-90,90,prm.num_users);
for i=1:prm.num_users
    userPos(1,i) = gNBpos(1,prm.refCellID) + r * cosd(theta(i));
    userPos(2,i) = gNBpos(2,prm.refCellID) + r * sind(theta(i));
end

userPos(1,1) = gNBpos(1,prm.refCellID) + r * cosd(10);
userPos(2,1) = gNBpos(2,prm.refCellID) + r * sind(10);

% Get scatters position
scatPos = get_scatters_positions(prm.numScat,prm.scen_center, ...
    prm.interSiteDist);

if prm.showFigures
    % Plot scenario
    show_scenario(gNBpos,prm.interSiteDist,cellCenters,prm.fillCell,userPos,scatPos)
end