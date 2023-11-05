function intPower = interf_gNBs_rx_power(prm,gNBpos,userPos,scatPos)

beamID = randi(prm.numTxBeams,prm.num_cells-1);
intPower = zeros(prm.num_users,prm.num_cells-1);
fprintf('<strong>Obtaining interfering RSRP of each user:</strong>\n');
for u=1:prm.num_users
    disp(['UserID: ' num2str(u)])
    for c=setdiff(1:prm.num_cells,prm.NCellID)
        % if c==12
        %     disp('cell 12')
        % end
        current_gNBpos = gNBpos(:,ceil(c/3));
        intPower(u,c) = interf_rxPowerPerUser(prm,current_gNBpos,c,userPos(:,u),...
            scatPos,beamID(c));
    end
end

min_db = -40;
process_db = @(x) max(min_db, x);

intPower = process_db(intPower);