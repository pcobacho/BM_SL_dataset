function [rxPower,RSRP,txBeamID,rxBeamID] = serving_gNB_rx_power(prm,gNBpos,userPos,scatPos)

RSRP = zeros(prm.numRxBeams,prm.numTxBeams,prm.num_users);
txBeamID = zeros(1,prm.num_users);
rxBeamID = zeros(1,prm.num_users);
rxPower = zeros(prm.num_users,1);
fprintf('<strong>Obtaining optimal beam pair RSRP of each user:</strong>\n');

current_gNBpos = gNBpos(:,ceil(prm.refCellID/3));

parfor u=1:prm.num_users
    disp(['UserID: ' num2str(u)])
        
    [RSRP(:,:,u),rxPower(u,1),txBeamID(1,u),rxBeamID(1,u)] = rxPowerPerUser(prm,current_gNBpos, ...
        prm.refCellID,userPos(:,u),scatPos);
end

if prm.showFigures
    IDs = sort(unique(txBeamID));
    colors = lines(length(IDs)+1);
    cIdx=1; % color index
    figure; hold on, grid on, axis equal
    for i=IDs
        uIdx = find(txBeamID==i); %user Index
        % plot(userPos(1,uIdx),userPos(2,uIdx),'Color',colors(cIdx,:),'Marker','o','LineStyle', 'none')
        scatter(userPos(1,uIdx),userPos(2,uIdx),30,colors(cIdx,:),'filled')
        cIdx = cIdx+1;
    end
end