% NR SSB Beam Sweeping
clear; clc; close all;
rng(211);   % Set RNG state for repeatability

% Add search path with the project files
% *************************************************************************
[~, oldPath] = addPaths();

%% SIMULATION PARAMETERS
prm.NCellID = 1;                    % Cell ID
prm.FreqRange = 'FR1';              % Frequency range: 'FR1' or 'FR2'
prm.CenterFreq = 3.2e9;              % Hz
prm.SSBlockPattern = 'Case C';      % Case A/B/C/D/E
prm.SSBTransmitted = [ones(1,8) zeros(1,0)];   % 4/8 or 64 in length

prm.TxArraySize = [8 8];            % Transmit array size, [rows cols]
prm.TxAZlim = [-180 -60];             % Transmit azimuthal sweep limits
prm.TxELlim = [-90 0];              % Transmit elevation sweep limits

prm.RxArraySize = [2 2];            % Receive array size, [rows cols]
prm.RxAZlim = [-180 180];           % Receive azimuthal sweep limits
prm.RxELlim = [0 90];               % Receive elevation sweep limits

prm.ElevationSweep = false;         % Enable/disable elevation sweep
prm.SNRdB = 30;                     % SNR, dB
prm.RSRPMode = 'SSSwDMRS';          % {'SSSwDMRS', 'SSSonly'}


prm = validateParams(prm);

%% SYNCHRONIZATION SIGNAL BURST CONFIGURATION
txBurst = nrWavegenSSBurstConfig;
txBurst.BlockPattern = prm.SSBlockPattern;
txBurst.TransmittedBlocks = prm.SSBTransmitted;
txBurst.Period = 20;
txBurst.SubcarrierSpacingCommon = prm.SubcarrierSpacingCommon;

% Configure an nrDLCarrierConfig object to use the synchronization signal
% burst parameters and to disable other channels. This object will be used
% by nrWaveformGenerator to generate the SS burst waveform.
cfgDL = configureWaveformGenerator(prm,txBurst);


%% BURST GENERATION
burstWaveform = nrWaveformGenerator(cfgDL);

% Display spectrogram of SS burst waveform
figure;
ofdmInfo = nrOFDMInfo(cfgDL.SCSCarriers{1}.NSizeGrid,prm.SCS);
nfft = ofdmInfo.Nfft;
spectrogram(burstWaveform,ones(nfft,1),0,nfft,'centered',ofdmInfo.SampleRate,'yaxis','MinThreshold',-130);
title('Spectrogram of SS burst waveform')


%% CHANNEL CONFIGURATION
c = physconst('LightSpeed');   % Propagation speed
lambda = c/prm.CenterFreq;     % Wavelength

prm.posTx = [0;0;0];           % Transmit array position, [x;y;z], meters
prm.posRx = [-30;-30;0];        % Receive array position, [x;y;z], meters

toRxRange = rangeangle(prm.posTx,prm.posRx);
spLoss = fspl(toRxRange,lambda);    % Free space path loss

% Transmit array (Conformal array for 360-degree sweep)  
numRows = prm.TxArraySize(1);
numCols = prm.TxArraySize(2);

C = arrayElementCoordinates(numRows,numCols,lambda);
theta = -120;
Cr = rotateArrayElements(C,theta,'Z');

taper = taylorwin(prm.NumTx);

arrayTx = phased.ConformalArray(...
    'ElementPosition', Cr', ...
    'Element', phased.IsotropicAntennaElement('BackBaffled', true),...
    'ElementNormal',[theta;0], ...
    'Taper',taper);

% figure, viewArray(arrayTx);
% 
% figure();
% pattern(arrayTx, prm.CenterFreq, -180:180, 0, 'PropagationSpeed', c,...
%     'CoordinateSystem', 'polar', ...
%     'Type', 'Directivity', 'PlotStyle', 'Overlay');

% Receive array
if prm.IsRxURA
    % Uniform rectangular array
    arrayRx = phased.URA(prm.RxArraySize,0.5*lambda, ...
        'Element',phased.IsotropicAntennaElement);
else
    % Uniform linear array
    arrayRx = phased.ULA(prm.NumRx, ...
        'ElementSpacing',0.5*lambda, ...
        'Element',phased.IsotropicAntennaElement);
end

% Scatterer locations
prm.enableScats = true;
prm.FixedScatMode = false;
if prm.enableScats
    if prm.FixedScatMode
        % Fixed single scatterer location
        prm.ScatPos = [50; 80; 0];     
    else
        % Generate scatterers at random positions
        Nscat = 20;        % Number of scatterers 
        azRange = -180:180;
        elRange = -90:90;
        radiusRange = 0:(abs(prm.posTx) + abs(prm.posRx));
    
        randAzOrder = randperm(length(azRange));
        randElOrder = randperm(length(elRange));
        ranRadOrder = randperm(length(radiusRange));
    
        azAngInSph = azRange(randAzOrder(1:Nscat));
        elAngInSph = elRange(randElOrder(1:Nscat));
        r = radiusRange(ranRadOrder(1:Nscat));    %radius
    
        [x,y,z] = sph2cart(deg2rad(azAngInSph),deg2rad(elAngInSph),r);
        prm.ScatPos = [x;y;z] + (prm.posTx + prm.posRx)/2;
    end
end

% Configure channel
channel = phased.ScatteringMIMOChannel;
channel.PropagationSpeed = c;
channel.CarrierFrequency = prm.CenterFreq;
channel.SampleRate = ofdmInfo.SampleRate;
channel.SimulateDirectPath = true;
channel.ChannelResponseOutputPort = true;
channel.Polarization = 'None';
channel.TransmitArray = arrayTx;
channel.TransmitArrayPosition = prm.posTx;
channel.ReceiveArray = arrayRx;
channel.ReceiveArrayPosition = prm.posRx;
channel.ScattererSpecificationSource = 'Property';
if isfield(prm,'ScatPos')
    channel.ScattererPosition = prm.ScatPos;
    channel.ScattererCoefficient = ones(1,size(prm.ScatPos,2));
end

% Get maximum channel delay
[~,~,tau] = channel(complex(randn(ofdmInfo.SampleRate*1e-3,prm.NumTx), ...
    randn(ofdmInfo.SampleRate*1e-3,prm.NumTx)));
maxChDelay = ceil(max(tau)*ofdmInfo.SampleRate);


%% TRANSMIT-END BEAM SWEEPING
% Number of beams at both transmit and receive ends
numBeams = sum(txBurst.TransmittedBlocks);

% Transmit beam angles in azimuth and elevation, equi-spaced
azBW = beamwidth(arrayTx,prm.CenterFreq,'Cut','Azimuth');
elBW = beamwidth(arrayTx,prm.CenterFreq,'Cut','Elevation');
txBeamAng = hGetBeamSweepAngles(numBeams,prm.TxAZlim,prm.TxELlim, ...
    azBW,elBW,prm.ElevationSweep);

% For evaluating transmit-side steering weights
SteerVecTx = phased.SteeringVector('SensorArray',arrayTx, ...
    'PropagationSpeed',c);

% Get the set of OFDM symbols occupied by each SSB
numBlocks = length(txBurst.TransmittedBlocks);
burstStartSymbols = ssBurstStartSymbols(txBurst.BlockPattern,numBlocks);
burstStartSymbols = burstStartSymbols(txBurst.TransmittedBlocks==1);
burstOccupiedSymbols = burstStartSymbols.' + (1:4);

% Apply steering per OFDM symbol for each SSB
gridSymLengths = repmat(ofdmInfo.SymbolLengths,1,cfgDL.NumSubframes);
%   repeat burst over numTx to prepare for steering
strTxWaveform = repmat(burstWaveform,1,prm.NumTx)./sqrt(prm.NumTx);
wT = zeros(prm.NumTx,numBeams);
for ssb = 1:numBeams

    % Extract SSB waveform from burst
    blockSymbols = burstOccupiedSymbols(ssb,:);
    startSSBInd = sum(gridSymLengths(1:blockSymbols(1)-1))+1;
    endSSBInd = sum(gridSymLengths(1:blockSymbols(4)));
    ssbWaveform = strTxWaveform(startSSBInd:endSSBInd,1);

    % Generate weights for steered direction
    wT(:,ssb) = SteerVecTx(prm.CenterFreq,txBeamAng(:,ssb));

    % Apply weights per transmit element to SSB
    strTxWaveform(startSSBInd:endSSBInd,:) = ssbWaveform.*(wT(:,ssb)');

end


%% RECEIVE-END BEAM SWEEPING AND MEASUREMENT
% Receive beam angles in azimuth and elevation, equi-spaced
azBW = beamwidth(arrayRx,prm.CenterFreq,'Cut','Azimuth');
elBW = beamwidth(arrayRx,prm.CenterFreq,'Cut','Elevation');
rxBeamAng = hGetBeamSweepAngles(numBeams,prm.RxAZlim,prm.RxELlim, ...
    azBW,elBW,prm.ElevationSweep);

% For evaluating receive-side steering weights
SteerVecRx = phased.SteeringVector('SensorArray',arrayRx, ...
    'PropagationSpeed',c);

% AWGN level
SNR = 10^(prm.SNRdB/20);                        % Convert to linear gain
N0 = 1/(sqrt(2.0*prm.NumRx*double(ofdmInfo.Nfft))*SNR); % Noise Std. Dev.

% Receive gain in linear terms, to compensate for the path loss
rxGain = 10^(spLoss/20);   

% Generate a reference grid for timing correction
%   assumes an SSB in first slot
carrier = nrCarrierConfig('NCellID',prm.NCellID);
carrier.NSizeGrid = cfgDL.SCSCarriers{1}.NSizeGrid;
carrier.SubcarrierSpacing = prm.SCS;
pssRef = nrPSS(carrier.NCellID);
pssInd = nrPSSIndices;
ibar_SSB = 0;
pbchdmrsRef = nrPBCHDMRS(carrier.NCellID,ibar_SSB);
pbchDMRSInd = nrPBCHDMRSIndices(carrier.NCellID);
pssGrid = zeros([240 4]);
pssGrid(pssInd) = pssRef;
pssGrid(pbchDMRSInd) = pbchdmrsRef;
refGrid = zeros([12*carrier.NSizeGrid ofdmInfo.SymbolsPerSlot]);
burstOccupiedSubcarriers = carrier.NSizeGrid*6 + (-119:120).';
refGrid(burstOccupiedSubcarriers, ...
    burstOccupiedSymbols(1,:)) = pssGrid;

% Loop over all receive beams
rsrp = zeros(numBeams,numBeams);
for rIdx = 1:numBeams

    % Fading channel, with path loss
    txWave = [strTxWaveform; zeros(maxChDelay,size(strTxWaveform,2))];
    fadWave = channel(txWave);

    % Receive gain, to compensate for the path loss
    fadWaveG = fadWave*rxGain;
    
    % Add WGN
    noise = N0*complex(randn(size(fadWaveG)),randn(size(fadWaveG)));
    rxWaveform = fadWaveG + noise;

    % Generate weights for steered direction
    wR = SteerVecRx(prm.CenterFreq,rxBeamAng(:,rIdx));

    % Apply weights per receive element
    if strcmp(prm.FreqRange, 'FR1')
        strRxWaveform = rxWaveform.*(wR');
    else  % for FR2, combine signal from antenna elements
        strRxWaveform = rxWaveform*conj(wR);
    end
    
    % Correct timing
    offset = nrTimingEstimate(carrier, ...
        strRxWaveform(1:ofdmInfo.SampleRate*1e-3,:),refGrid*wR(1)');
    if offset > maxChDelay
        offset = 0;
    end
    strRxWaveformS = strRxWaveform(1+offset:end,:);

    % OFDM Demodulate
    rxGrid = nrOFDMDemodulate(carrier,strRxWaveformS);

    % Loop over all SSBs in rxGrid (transmit end)
    for tIdx = 1:numBeams
        % Get each SSB grid
        rxSSBGrid = rxGrid(burstOccupiedSubcarriers, ...
            burstOccupiedSymbols(tIdx,:),:);

        if strcmpi(prm.RSRPMode,'SSSwDMRS')
            meas = nrSSBMeasurements(rxSSBGrid,carrier.NCellID,mod(tIdx-1,8));
        else
            meas = nrSSBMeasurements(rxSSBGrid,carrier.NCellID);
        end
        rsrp(rIdx,tIdx) = max(meas.RSRPPerAntenna);
    end
end


%% BEAM DETERMINATION
[m,i] = max(rsrp,[],'all','linear');    % First occurrence is output
% i is column-down first (for receive), then across columns (for transmit)
[rxBeamID,txBeamID] = ind2sub([numBeams numBeams],i(1));

% Display the selected beam pair
disp(['Selected Beam pair with RSRP: ' num2str(rsrp(rxBeamID, ...
    txBeamID)) ' dBm', 13 '  Transmit #' num2str(txBeamID) ...
    ' (Azimuth: ' num2str(txBeamAng(1,txBeamID)) ', Elevation: ' ...
    num2str(txBeamAng(2,txBeamID)) ')' 13 '  Receive #' num2str(rxBeamID) ...
    ' (Azimuth: ' num2str(rxBeamAng(1,rxBeamID)) ', Elevation: ' ...
    num2str(rxBeamAng(2,rxBeamID)) ')' ]);

% Display final beam pair patterns
h = figure('Position',figposition([32 55 32 40]));
h.Name = 'Selected Transmit Array Response Pattern';
wT_opt = SteerVecTx(prm.CenterFreq,txBeamAng(:,txBeamID));
pattern(arrayTx,prm.CenterFreq,'PropagationSpeed',c,'Weights',wT_opt);

h = figure('Position',figposition([32 55 32 40]));
h.Name = 'Selected Receive Array Response Pattern';
wR = SteerVecRx(prm.CenterFreq,rxBeamAng(:,rxBeamID));
pattern(arrayRx,prm.CenterFreq,'PropagationSpeed',c,'Weights',wR);

h = figure();
subplot(1,2,1)
pattern(arrayTx, prm.CenterFreq, -180:180, 0, 'PropagationSpeed', c,...
'CoordinateSystem', 'polar', 'Type', 'Directivity', 'PlotStyle', 'Overlay',...
'Weights',wT(:,txBeamID))

subplot(1,2,2)
pattern(arrayRx, prm.CenterFreq, -180:180, 0, 'PropagationSpeed', c,...
'CoordinateSystem', 'polar', 'Type', 'Directivity', 'PlotStyle', 'Overlay',...
'Weights',wR)

% Plot MIMO scenario with tx, rx, scatterers, and determined beams. Beam
% patterns in this figure resemble the power patterns in linear scale.
prmScene = struct();
prmScene.TxArray = arrayTx;
prmScene.RxArray = arrayRx;
prmScene.TxArrayPos = prm.posTx;
prmScene.RxArrayPos = prm.posRx;
if isfield(prm,'ScatPos')
    prmScene.ScatterersPos = prm.ScatPos;
end
prmScene.Lambda = lambda;
prmScene.ArrayScaling = 1;     % To enlarge antenna arrays in the plot
prmScene.MaxTxBeamLength = 32; % Maximum length of transmit beams in the plot
prmScene.MaxRxBeamLength = 10; % Maximum length of receive beam in the plot
hPlotSpatialMIMOScene(prmScene,wT(:,txBeamID),wR);
if ~prm.ElevationSweep
    view(2);
end

% Plot the scattering MIMO scenario (including transmit and receive antenna
% arrays, scatterer positions and their paths, and all the transmit and 
% receive antenna array beam patterns) by using the helper function 
% hPlotSpatialMIMOScene. Beam patterns in this figure resemble the power
% patterns in linear scale.Plot the scattering MIMO scenario (including
% transmit and receive antenna arrays, scatterer positions and their paths,
% and all the transmit and receive antenna array beam patterns) by using 
% the helper function hPlotSpatialMIMOScene. Beam patterns in this figure
% resemble the power patterns in linear scale.Plot the scattering MIMO
% scenario (including transmit and receive antenna arrays, scatterer 
% positions and their paths, and all the transmit and receive antenna array
% beam patterns) by using the helper function hPlotSpatialMIMOScene. Beam 
% patterns in this figure resemble the power patterns in linear scale.
hPlotSpatialMIMOScene(prmScene,wT,wR);
axis tight;
view([74 29]);

% Restore search paths
% *************************************************************************
path(oldPath);