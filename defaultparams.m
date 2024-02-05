function prm = defaultparams()

% Scenario parameters
prm.scen_center = [0,0];  % center of scenario
prm.num_users = 20;      % number of UEs in reference cell
prm.hUE = 1.5;              % UEs height
prm.num_gNB = 7;          % number of base stations (gNBs)
prm.h_gNB = 25;            % gNB height
prm.num_cells = prm.num_gNB*3; % number of cells (sectors)
prm.interSiteDist = 200;  % inter site distance
prm.refCellID = 1;        % reference cell ID (cell under study)
prm.numScat = 0;         % number of scatters in scenario
prm.fillCell = false;    % fills cells with color according to their gNB

% Configure Antenna Array, Frequency and Beam Sweeping Angles
prm.NCellID = 1;               % Cell ID
prm.FreqRange = 'FR1';              % Frequency range: 'FR1' or 'FR2'
prm.CenterFreq = 3.2e9;              % Hz
prm.SSBlockPattern = 'Case C';      % Case A/B/C/D/E
prm.SSBTransmitted = [ones(1,8) zeros(1,0)];   % 4/8 or 64 in length

prm.TxArraySize = [8 8];            % Transmit array size, [rows cols]
prm.TxAZranges = [[-50 50]; [70 170]; [-170 -70]];
prm.TxELlim = [-60 -15];              % Transmit elevation sweep limits

prm.RxArraySize = [2 2];            % Receive array size, [rows cols]
prm.RxAZlim = [0 180];           % Receive azimuthal sweep limits
prm.RxELlim = [0 90];               % Receive elevation sweep limits
prm.rxGain_dB = 98;                 % Rx antenna gain (in dB)

prm.cellAngles = [0, 120, -120];

prm.ElevationSweep = false;         % Enable/disable elevation sweep
prm.SNRdB = 30;                     % SNR, dB
prm.RSRPMode = 'SSSwDMRS';          % {'SSSwDMRS', 'SSSonly'}

% Simulation parameters
prm.showFigures = false;
prm.resultFolder = 'sim';
prm.saveResults = true;