function intPower = interf_gNBs_rx_power(prm,gNBpos,userPos,scatPos)

beamID = randi([1,prm.numTxBeams],1,prm.num_cells);
intPower = zeros(prm.num_users,prm.num_cells-1);
cells = setdiff(1:prm.num_cells,prm.NCellID);

fprintf('<strong>Obtaining interfering RSRP of each user:</strong>\n');
parfor u=1:prm.num_users
    disp(['UserID: ' num2str(u)])
    for c=cells
        current_gNBpos = gNBpos(:,ceil(c/3));
        intPower(u,c) = interf_rxPowerPerUser(prm,current_gNBpos,c,userPos(:,u),...
            scatPos,beamID(c));
    end
end

min_db = -40;
process_db = @(x) max(min_db, x);

intPower = process_db(intPower);