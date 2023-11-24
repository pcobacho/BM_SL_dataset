% *************************************************************************
% AUTHORS: Francisco J. Martin Vega
% *************************************************************************
% GROUP: Lab 1.3.5., Communications and Signal Processing Lab (ComSP), 
% TELMA, ETSIT, University of Malaga
% *************************************************************************
% DESCRIPTION:
% This script creates the results folder
% *************************************************************************

function filename = createFilename(prm,seed)

numTxBeams = length(prm.SSBTransmitted);
path = [cd '\Results'];
filename = sprintf('%dusers_%dbeams_%dscat_seed%d.mat',prm.num_users,numTxBeams,prm.numScat,seed);

full_path = fullfile(path,filename);

flnameAux = full_path;
cont=1;
while isfile(flnameAux)    
    aux = sprintf('%dusers_%dbeams_%dscat_seed%d(%d).mat',prm.num_users,...
        numTxBeams,prm.numScat,seed,cont);
    flnameAux = fullfile(path,aux);
    cont=cont+1;
end

filename = flnameAux;