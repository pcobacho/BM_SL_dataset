function [rsrp,optRSRP,txBeamID,rxBeamID] = rxPowerPerUser(prm,gNBpos,userPos,scatPos)

%% SIMULATION PARAMETERS
TxAZlim = prm.TxAZranges(mod(prm.NCellID-1,3)+1,:);  % Transmit azimuthal sweep limits

theta = prm.cellAngles(mod(prm.NCellID-1,3)+1); % Current cell normal angle

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
% figure;
ofdmInfo = nrOFDMInfo(cfgDL.SCSCarriers{1}.NSizeGrid,prm.SCS);
nfft = ofdmInfo.Nfft;

%% CHANNEL CONFIGURATION
c = physconst('LightSpeed');   % Propagation speed
lambda = c/prm.CenterFreq;     % Wavelength

toRxRange = rangeangle(gNBpos,userPos);
spLoss = fspl(toRxRange,lambda);    % Free space path loss

% Transmit array (Conformal array for 360-degree sweep)  
numRows = prm.TxArraySize(1);
numCols = prm.TxArraySize(2);

C = arrayElementCoordinates(numRows,numCols,lambda);
Cr = rotateArrayElements(C,theta,'Z');

taper = taylorwin(prm.NumTx);

arrayTx = phased.ConformalArray(...
    'ElementPosition', Cr', ...
    'Element', phased.IsotropicAntennaElement('BackBaffled', true),...
    'ElementNormal',[theta;0], ...
    'Taper',taper);

arrayRx = phased.ConformalArray('Element', ...
    phased.IsotropicAntennaElement('BackBaffled', false));

% Scatterer locations
if ~isempty(scatPos)
    prm.ScatPos = scatPos;
end

% Configure channel
switch prm.chType
    case 'scat'
        channel = phased.ScatteringMIMOChannel;
        channel.PropagationSpeed = c;
        channel.CarrierFrequency = prm.CenterFreq;
        channel.SampleRate = ofdmInfo.SampleRate;
        channel.SimulateDirectPath = true;
        channel.ChannelResponseOutputPort = true;
        channel.Polarization = 'None';
        channel.TransmitArray = arrayTx;
        channel.TransmitArrayPosition = gNBpos;
        channel.ReceiveArray = arrayRx;
        channel.ReceiveArrayPosition = userPos;
        channel.ScattererSpecificationSource = 'Auto';
        channel.NumScatterers = prm.numScat;
        if isfield(prm,'ScatPos')
            channel.ScattererPosition = prm.ScatPos;
            channel.ScattererCoefficient = ones(1,size(prm.ScatPos,2));
        end
        
        % Get maximum channel delay
        [~,~,tau] = channel(complex(randn(ofdmInfo.SampleRate*1e-3,prm.NumTx), ...
            randn(ofdmInfo.SampleRate*1e-3,prm.NumTx)));
        maxChDelay = ceil(max(tau)*ofdmInfo.SampleRate);
    case 'freeSp'
        channel = phased.FreeSpace('SampleRate',ofdmInfo.SampleRate, ...
            'PropagationSpeed',c,'OperatingFrequency',prm.CenterFreq);
        maxChDelay = 0;
    otherwise
        error('prm.chType is invalid. Must be ''scat'' or ''freeSp''.')
end


%% TRANSMIT-END BEAM SWEEPING

if prm.numTxBeams>1
    % Transmit beam angles in azimuth and elevation, equi-spaced
    % azBW = beamwidth(arrayTx,prm.CenterFreq,'Cut','Azimuth');
    % elBW = beamwidth(arrayTx,prm.CenterFreq,'Cut','Elevation');
    % txBeamAng = hGetBeamSweepAngles(prm.numTxBeams,TxAZlim,prm.TxELlim, ...
    %     azBW,elBW,prm.ElevationSweep);

    txBeamAng = [linspace(TxAZlim(1),TxAZlim(2),prm.numTxBeams); zeros(1,prm.numTxBeams)];
else
    txBeamAng = [0;0];
end

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
% repeat burst over numTx to prepare for steering
strTxWaveform = repmat(burstWaveform,1,prm.NumTx)./sqrt(prm.NumTx);
wT = zeros(prm.NumTx,prm.numTxBeams);

for ssb = 1:prm.numTxBeams
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

reset(SteerVecTx)
release(SteerVecTx)

%% RECEIVE-END BEAM SWEEPING AND MEASUREMENT

% Receive beam angles in azimuth and elevation, equi-spaced
azBW = beamwidth(arrayRx,prm.CenterFreq,'Cut','Azimuth');
elBW = beamwidth(arrayRx,prm.CenterFreq,'Cut','Elevation');
rxBeamAng = hGetBeamSweepAngles(prm.numRxBeams,prm.RxAZlim,prm.RxELlim, ...
    azBW,elBW,prm.ElevationSweep);

% For evaluating receive-side steering weights
SteerVecRx = phased.SteeringVector('SensorArray',arrayRx, ...
    'PropagationSpeed',c);

% AWGN level
SNR = 10^(prm.SNRdB/20);                        % Convert to linear gain
N0 = 1/(sqrt(2.0*prm.NumRx*double(ofdmInfo.Nfft))*SNR); % Noise Std. Dev.

% Receive gain in linear terms, to compensate for the path loss
% rxGain = 10^(spLoss/20);   
rxGain = 10^(prm.rxGain_dB/20);

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
rsrp = zeros(prm.numRxBeams,prm.numTxBeams);
for rIdx = 1:prm.numRxBeams
    txWave = [strTxWaveform; zeros(maxChDelay,size(strTxWaveform,2))];

    switch prm.chType
        case 'scat'
            fadWave = channel(txWave);
        case 'freeSp'
            fadWave = channel(txWave,Cr.',userPos,zeros(3,prod(prm.TxArraySize)),[0;0;0]);
            fadWave = sum(fadWave,2);
    end

    % Receive gain, to compensate for the path loss
    fadWaveG = fadWave*rxGain;
    
    % no added AWGN
    % noise = N0*complex(randn(size(fadWaveG)),randn(size(fadWaveG)));
    rxWaveform = fadWaveG;% + noise;

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
    for tIdx = 1:prm.numTxBeams
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
[~,i] = max(rsrp,[],'all','linear');    % First occurrence is output
% i is column-down first (for receive), then across columns (for transmit)
[rxBeamID,txBeamID] = ind2sub([prm.numRxBeams prm.numTxBeams],i(1));
optRSRP = rsrp(rxBeamID,txBeamID);