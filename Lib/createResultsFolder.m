% *************************************************************************
% AUTHORS: Francisco J. Martin Vega
% *************************************************************************
% GROUP: Lab 1.3.5., Communications and Signal Processing Lab (ComSP), 
% TELMA, ETSIT, University of Malaga
% *************************************************************************
% DESCRIPTION:
% This script creates the results folder
% *************************************************************************

function prm = createResultsFolder(prm)

[~, message, ~] = mkdir('Results');

strDate = char(datetime('today'));
prm.folderName = ['Results/' prm.resultFolder '_' strDate];
[~, message, ~] = mkdir(prm.folderName);
i = 1;
folderNameAux = prm.folderName;
while strcmp(message, 'Directory already exists.')
    folderNameAux = [prm.folderName ' (' num2str(i) ')'];
    i = i + 1;
    [~, message, ~] = mkdir(folderNameAux);
end
prm.folderName = folderNameAux;

