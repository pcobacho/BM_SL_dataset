% *************************************************************************
% AUTHORS: Pablo Cobacho
% *************************************************************************
% GROUP: Lab 1.3.5., Communications and Signal Processing Lab (ComSP), 
% TELMA, ETSIT, University of Malaga
% *************************************************************************
% DESCRIPTION:
% This script simulates the P2 beam management procedure in a mresobile 
% communications scenario with hexagonal cells. In addition, it calculates
% the downlink SINR for the users contained in the cell under study 
% (refCellID).
% *************************************************************************

clear; clc; close all;
rng(211);   % Set RNG state for repeatability
tic;
% Add search path with the project files
[~, oldPath] = addPaths();

% Default parameters
prm = defaultparams();

% Customize parameters
prm.num_users = 1;
prm.fillCell = false;
prm.ElevationSweep = false;

prm.numScat = 0;
prm.interSiteDist = 200;

prm.hUE = 0;              % UEs height
prm.h_gNB = 0;

% prm.FreqRange = 'FR2'; % ¿No funciona?
% prm.CenterFreq = 32e9;
% prm.SSBlockPattern = 'Case D';
% prm.SSBTransmitted = [ones(1,64) zeros(1,0)];

% Parameters validation
prm = validateParams(prm);

% Scenario initialization
[gNBpos,cellCenters,userPos,scatPos] = iniScenario(prm);

% Obtains the optimal beam pair RSRP each user of the reference cell
rxPowerdBm = serving_gNB_rx_power(prm,gNBpos,userPos,scatPos);
rxPowerdB = rxPowerdBm - 30;

% Calculates interfering RSRP
intPowerdBm = interf_gNBs_rx_power(prm,gNBpos,userPos,scatPos);
% intPowerdBm = intPowerdBm(:,2:end);
intPowerdB = intPowerdBm-30;
intPower = 10.^(intPowerdB./20);
totalIntPower = sum(intPower,2);
totalIntPowerdB = 20*log10(totalIntPower);

% Calculate SINR per user
BWinHz = 52*12*30e3;
F=3; % noise figure
NdB = -173.8+10*log10(BWinHz)-30+F;

rxPower =  10.^(rxPowerdB./20);
N = 10^(NdB/20);

sinrdB = 10*log10(rxPower./(totalIntPower+N));

% Restore search paths
% *************************************************************************
path(oldPath);
toc;