function rxPower = serving_gNB_rx_power(prm,gNBpos,userPos,scatPos)

txBeamID = zeros(1,prm.num_users);
rxBeamID = zeros(1,prm.num_users);
rxPower = zeros(prm.num_users,1);
fprintf('<strong>Obtaining optimal beam pair RSRP of each user:</strong>\n');
for u=1:prm.num_users
    disp(['UserID: ' num2str(u)])
    
    current_gNBpos = gNBpos(:,ceil(prm.refCellID/3));
    [~,rxPower(u,1),txBeamID(1,u),rxBeamID(1,u)] = rxPowerPerUser(prm,current_gNBpos, ...
        prm.refCellID,userPos(:,u),scatPos);
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