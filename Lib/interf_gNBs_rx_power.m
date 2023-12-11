function intPower = interf_gNBs_rx_power(prm,gNBpos,userPos,scatPos,rxBeamID)

beamID = randi([1,prm.numTxBeams],1,prm.num_cells);
intPower = zeros(prm.num_users,prm.num_cells-1);
cells = setdiff(1:prm.num_cells,prm.NCellID);

fprintf('<strong>Obtaining interfering RSRP of each user:</strong>\n');
for c=cells
    fprintf('CellID: %d\n', c);
    current_gNBpos = gNBpos(:,ceil(c/3));
    bID = beamID(c);
    for u=1:prm.num_users
        fprintf('\tUserID: %d\n', u);
        intPower(u,c) = interf_rxPowerPerUser(prm,current_gNBpos,c,userPos(:,u),...
            scatPos,bID,rxBeamID(u));
    end
end

min_db = -40;
process_db = @(x) max(min_db, x);

intPower = process_db(intPower);