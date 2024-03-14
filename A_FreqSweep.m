%% Sweep in laser frequency
% Might change fcf, needed for multiring structures

% Structure that's being worked on
structure = "newChip_ring1";

% Splitter at input, value
splitter = 99/1; %50/50;  % 99/1

% Detector self-amplifications
boostPDA = 10^(20/20);
boostWX = 50;

% Choose detector calibration
%detectorFolder = 'Instrument_Characterization\Instrument_data\CalibrateDetectors\DetectorData2023-06-16_16-32.mat';
% detectorFolder = 'Instrument_Characterization\Instrument_data\CalibrateDetectors\DetectorData2023-07-31_10-42.mat';

% Comments that are added to final savefile
conditions = "Cleanmode(0)";

%% Convert values chosen above to measurement parameters
% cd 'C:\Users\Nanolab\Desktop\DavideM\'
% addpath 'Functions'
saveFolder = 'Saved_data/A_sweepFreq/';
[~,~] = mkdir(saveFolder);

load('DetectorData2023-07-31_10-42.mat')
gainWX = gainWX*boostWX;
gainPDA = gainPDA*boostPDA;

% clear daqCurrents
switch structure
    case "R-1"
        f0 = 192.065;  %Thz
        %freq_range = f0 + (-30e-3:0.5e-3:30e-3);
        freq_range = f0 + (-30e-3:1e-3:30e-3);

    case "R-2"
        %f0 = 192.663;   %THz
            % Old setup with Peltier
        f0 = 192.653;
        freq_range = f0 + (-30e-3:0.5e-3:30e-3);
        
    case "R-3"
        f0 = 192.655;   %THz
        freq_range = f0 + (-30e-3:0.5e-3:30e-3);
    case "CROW-2-1"
        freq_range = (193.536000 : 0.5e-3 : 193.67500);

        daqCurrents.addoutput("DAC_Current", 0, "Current");
            % These are ring bot mid and top, for now 
        daqCurrents.write(0)
    case "CROW-3-2"
        f0 = 192.59000;
        freq_range = f0 + (-30e-3:0.5e-3:150e-3);
    case "CROW-2-1-heated"
        f0 = 191.77000;
        freq_range = f0 + (-30e-3:0.5e-3:300e-3);

        daqCurrents = instrumentation.heaterDaq.heaterDaq(0:2);
            % These are ring bot mid and top, for now 
        daqCurrents.currents([18,0,0]);
    case "SCISSOR-2-1"
        freq_range = (192.02000 : 1e-3 : 192.115);

        daqCurrents.addoutput("DAC_Current", 1:2, "Current");
        daqCurrents.write([0,0])
    case "SCISSOR-2-1_14mA"
        freq_range = (192.02000 : 1e-3 : 192.115);

        
        daqCurrents = instrumentation.heaterDaq.heaterDaq(1:2);
            % These are ring bot mid and top, for now 
        daqCurrents.currents([14,0]);
    case "SCISSOR-2-1_16mA"
        freq_range = (192.02000 : 1e-3 : 192.09);

        daqCurrents = instrumentation.heaterDaq.heaterDaq(1:2);
            % These are ring bot mid and top, for now 
        daqCurrents.currents([16,0]);
    case "SCISSOR-3-2"
        freq_range = (191.5 : 3e-3 : 192.5);
    case "newChip_SCISSOR-2-1"
        freq_range = 192.360 : 1e-3: 192.47;
    case "newChip_SCISSOR-2-1_10mA"
        % freq_range = 192.355 : 1e-3: 192.450;
        freq_range = 192.1 : 1e-2: 193.450;
    case "newChip_SCISSOR-3-2"
        freq_range = 192.38 : 1e-3: 192.46;
        % freq_range = 192.20 : 1e-2: 192.7;
        % freq_range = 192.0 : 1e-1: 194;
    case "newChip_SCISSOR-1-3"
        % freq_range = 192.0 : 1e-2: 194.0;   % resonance around 193 THz
        freq_range = 192.96 : 1e-3: 193.04;
        % daqCurrents = instrumentation.heaterDaq.heaterDaq(1:2);
            % These are ring bot mid and top, for now 
        % daqCurrents.currents([10,0]);
    case "newChip_SCISSOR-3-3"
        % freq_range = 192.30 : 1e-2: 192.70;
        freq_range = 192.37 : 1e-2: 192.47;
    case "newChip_SCISSOR-3-MultiOut"
        % freq_range = 192.35 : 0.5e-2: 192.55;
        freq_range = 192.41 : 1e-3: 192.51;
    case "newChip_ring1"
        freq_range = 192.30 : 2e-3: 192.46;
        
end

%% Initialize instruments
% picoON;
pico.set_closest_sampling_time( 1e-4 );
pico.meas_time(0.5);
pico.SetChannel(enable = [1,1,1,1]);

% VOAin_on;
% VOAout_on;

f0 = min(freq_range)+30e-3;
% tlsON;


% Initialize output
drop = nan(length(freq_range),pico.samples_per_trace);
through = nan(length(freq_range),pico.samples_per_trace);
input = nan(length(freq_range),1);


%% Calculate zeroPoints of detectors

TLS.disable;
VOAin(0);
VOAout(0);
TLS.wait_pending_operation;

pause(0.5)
pico.acquire(autoRange = 1);

zeroFast = mean(pico.channels.data{1});
zeroPDA = mean(pico.channels.data{2});
zeroWX = mean(pico.channels.data{3});
disp([zeroFast, zeroPDA, zeroWX]);


%% Perform the measurement
% Liveplot
fig = figure(Name = "Frequency sweep", Units = "normalized", Position=[0.1,0.1,0.5,0.5]);
axs = axes(Parent=fig, Units='normalized');
hold(axs,'on');
grid(axs, 'on'); grid(axs, 'minor');
axs.XLabel.String = 'Frequency (THz)';
axs.YLabel.String = 'Avg power (au)';

% Maybe add errors?
pl(1) = plot(axs, freq_range, mean(drop,2), color="blue", DisplayName = "Drop",...
    LineStyle="-", Marker=".", YDataSource="(mean(drop,2)-zeroFast)/gainFast./ ((input-zeroWX)/gainWX*splitter)");
pl(2) = plot(axs, freq_range, mean(through,2), color="red", DisplayName = "Through", ...
    LineStyle="-", Marker=".", YDataSource="(mean(through,2)-zeroPDA)/gainPDA./ ((input-zeroWX)/gainWX*splitter)");
legend(location="east");
refreshdata;
xlim("tight")

% axs.XLim = [min(pl(1).XData), max(pl(1).XData)];

TLS.cleanMode(0);
TLS.ftf(-30e3);
TLS.enable;

VOAin(0.2);
VOAout(0.02);

for i=1:length(freq_range)
    % Get things right converts to unsigned int for the LASER
    
    TLS.freq(freq_range(i));
    refreshdata; drawnow;

    TLS.wait_pending_operation;
    pause(0.1);

    pico.acquire(autoRange = 1);
    % AdjustFastDetector;
    % Output to channel B

    pico.acquire(autoRange=true);

    drop(i,:) = (pico.channels.data{4}-zeroFast)+zeroFast;
    through(i,:) = pico.channels.data{1};
    input(i) = mean(pico.channels.data{3});


    if(max(through)>3.5)
        VOAin.pwr(0.0);
        error("saturating")
    end
    disp(i+"/"+length(freq_range));
end
VOAin(0.02);
VOAout(0.1);

timestamp = datetime; % Save a timestamp, makes retrocompatibily handling much easier
save(saveFolder + "Structure_"+structure+"_"+ string(datetime("now","format",'yyyy-MM-dd_HH-mm')),...
    "conditions", "freq_range", "drop", "through", "input", "gainFast",...
    "gainWX", "gainPDA", "zeroFast", "zeroWX", "zeroPDA", "splitter",...
    "timestamp");


%% Plot result
% load(detectorFolder)

% Legacy
if(~exist("freq_range","var"))
    freq_range = f0 + det_range/10^6;
end
if(exist("timestamp", "var"))
    if(timestamp<datetime("16-Jun-2023 14:35:10"))
        
    end
else
    % Correct the coGain (wrong coupling on WX)
    gainWX = gainWX*1.577;
    gainFast = gainFast*1.577;
    gainPDA = gainPDA*1.577;

    boostPDA = 1;
    boostWX = nan;
end

inputW = (input-zeroWX)/gainWX*splitter;
throughW = (mean(through,2)-zeroPDA)/gainPDA;
dropW = (mean(drop,2)-zeroFast)/gainFast;

% % % OVERRIDE
%  inputW = 1;
%  throughW = mean(through,2);
%  dropW = mean(drop,2);





fig = figure(Name = "Frequency sweep", Position = [100,100,1200,400]);
axs = axes(Parent=fig, Units='normalized');
hold(axs,'on');
grid(axs, 'on'); grid(axs, 'minor');
axs.XLabel.String = 'Frequency (Thz)';
axs.YLabel.String = 'Transmission';

% First plot
x = transpose(freq_range); % Frequency in THz
y_th = throughW./inputW;
y_dr = dropW./inputW;




pl(1) = plot(axs, x, y_th, color="red", DisplayName = "Through",...
    LineStyle="-", Marker=".");
pl(2) = plot(axs, x, y_dr, color="blue", DisplayName = "Drop", ...
    LineStyle="-", Marker=".");

axs.XLim = [min(pl(1).XData), max(pl(1).XData)];
%ylim([-0.0002, max(y_th)+0.0002]);

leg = legend(location="east");

VOAout(0.1);   % lower this to be sure not to overpower the PD on the through port

% %% Microrings: add a lorentian fit and plot it
% if(structure == "R-1" || structure == "R-2" || structure == "R-3")
%     % For microrings, fit lorentian
%     peakData_through = fitLorentian(x,y_th);
%     peakData_drop = fitLorentian(x,y_dr);
%     pl(3) = plot(axs, x, peakData_drop(x),...
%     color="green", DisplayName = "Best fit drop", linestyle = "--");
%     pl(4) = plot(axs, x, peakData_through(x),...
%     color="magenta", DisplayName = "Best fit through", linestyle = "--");
% 
%     disp("Measured quality factors = " + peakData_through.ctr/peakData_through.HWHM + "  " + peakData_drop.ctr/peakData_drop.HWHM)
%     fprintf("Measured center frequencies = %6.7f   %6.7f", peakData_through.ctr, peakData_drop.ctr)
%     disp(" ")
% end

