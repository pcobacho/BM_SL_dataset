% *************************************************************************
% AUTHORS: Pablo Cobacho
% *************************************************************************
% GROUP: Lab 1.3.5., Communications and Signal Processing Lab (ComSP), 
% TELMA, ETSIT, University of Malaga
% *************************************************************************
% DESCRIPTION:
% This script simulates the P2 beam management procedure in a mobile 
% communications scenario with hexagonal cells. In addition, it calculates
% the downlink SINR for the users contained in the cell under study 
% (refCellID).
% *************************************************************************

clear; clc; close all;
seed=219;
rng(seed);   % Set RNG state for repeatability
tic;

% Add search path with the project files
[~, oldPath] = addPaths();

% Default parameters
prm = defaultparams();

% Customize parameters
prm.num_users = 13;
prm.fillCell = false;

prm.numScat = 0;
prm.interSiteDist = 200;

prm.hUE = 0;              % UEs height
prm.h_gNB = 0;
prm.ElevationSweep = false;

prm.FreqRange = 'FR2';
prm.CenterFreq = 32e9;
prm.SSBlockPattern = 'Case D';
numTxBeams = 1;
prm.SSBTransmitted = [ones(1,numTxBeams) zeros(1,64-numTxBeams)];

prm.showFigures = true;
prm.saveResults = true;

% Parameters validation
prm = validateParams(prm);

% Scenario initialization
[gNBpos,cellCenters,userPos,scatPos] = iniScenario(prm);

% Obtains the optimal beam pair RSRP each user of the reference cell
[rxPowerdBm,RSRP,txBeamID,rxBeamID] = serving_gNB_rx_power(prm,gNBpos,userPos,scatPos);
rxPowerdB = rxPowerdBm - 30;

% Calculate SINR per user
BWinHz = 52*12*prm.SCS*1e3;
F=3; % noise figure
NdB = -173.8+10*log10(BWinHz)-30+F;

rxPower =  10.^(rxPowerdB./10);
N = 10^(NdB/10);

% sinrdB = 10*log10(rxPower./(totalIntPower+N));
sinrdB = 10*log10(rxPower./N);

%% BORRAR
rsrp = squeeze(RSRP);
angle = linspace(-90,90,90);
figure, polarplot(deg2rad(angle),rsrp.','.-')
%%

if prm.showFigures
    % SINR histogram (pdf)
    figure, histogram(sinrdB,'Normalization','pdf')
    xlabel('SINR (dB)'), title('SINR PDF')

    % SINR map
    if prm.num_users>1
        showUserSINRmap(sinrdB,userPos,prm.num_users)
    end
end

if prm.saveResults    
    % Save results
    filename = createFilename(prm,seed);
    save(filename,'prm','sinrdB','RSRP','txBeamID','rxBeamID');
end

% Restore search paths
% *************************************************************************
path(oldPath);

toc;