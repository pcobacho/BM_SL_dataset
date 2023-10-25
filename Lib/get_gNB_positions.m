% *************************************************************************
% AUTHORS: Pablo Cobacho
% *************************************************************************
% GROUP: Lab 1.3.5., Communications and Signal Processing Lab (ComSP), 
% TELMA, ETSIT, University of Malaga
% *************************************************************************
% DESCRIPTION:
% this script calculates the coordinates of each of the base stations in
% the scenario. The first cell is placed in the centre of the stage 
% (scen_center) and the other cells are placed around the first cell in
% a ring shape.
% *************************************************************************

function  gNBpos = get_gNB_positions(scen_center,interSiteDist,gNB_height)

num_gNB = 7;  % number of gNBs
gNB_ang = linspace(0, 360, num_gNB);  % gNB position angles (in degrees)
gNB_ang = gNB_ang(1:end-1);

gNBpos(1,:) = [scen_center(1), scen_center(1) + interSiteDist * cosd(gNB_ang)];
gNBpos(2,:) = [scen_center(2), scen_center(2) + interSiteDist * sind(gNB_ang)];
gNBpos(3,:) = gNB_height*ones(1,num_gNB);
