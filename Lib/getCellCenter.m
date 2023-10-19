% *************************************************************************
% AUTHORS: Pablo Cobacho
% *************************************************************************
% GROUP: Lab 1.3.5., Communications and Signal Processing Lab (ComSP), 
% TELMA, ETSIT, University of Malaga
% *************************************************************************
% DESCRIPTION: 
% From the positions of the base stations in the scenario (gNBpos) and the 
% spacing between them (interSiteDist), this script calculates the 
% coordinates of the three cell centres in each trisector base station.
% *************************************************************************

function cellCenters = getCellCenter(gNBpos,interSiteDist)

num_gNB = length(gNBpos);
cellCenters = zeros(2,3,num_gNB);

r = interSiteDist/3; % Cell radius
cell_angles = [0 120 240];

for i=1:num_gNB
    cellCenters(1,:,i) = gNBpos(1,i) + r * cosd(cell_angles);
    cellCenters(2,:,i) = gNBpos(2,i) + r * sind(cell_angles);
end