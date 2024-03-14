%% Here is a script to the self pulsing map
% transmission at minimum value
% for 3-rings scissor 1 (close to the 2-ring scissor): VOAin pwr0.1 gives more
% or less linear, 0.2 is already crazy only for freq 192.43. For pwr0.3.
% self pulsing from 192.40 to 192.44. At 192.42, SP only from pwr 0.4.
% --> good range is pwr from 0.1 to 0.6, freq from 192.40 to 192.44
% NB: set oscilloscope with coupling range (input impedance) of 50 ohm!

% Structure that's being worked on
structure = "structure_1";
cd 'C:\Users\Nanolab\Desktop\Alessio_Lugnan\Chip_PELM_DavideMicheli\SCISSOR_2rings\Matlab_functions_and_scripts\'
saveFolder = 'Saved_data/1R_SPmap_class/';

% initalize instruments
% Init_VOAKeySight_in_out
% tlsON;
% awgON;
% picoON;

% parameters
gain_PDin = 2e3;
EDFA_current = 400.9   % (mA)
saturation_voltage = 1.6;   % saturation voltage level of the sensitive PD, usually at the drop port

switch structure
    case "newChip_SCISSOR-3-1"
        laserfreq_sweep = 192.40:0.005:192.44;   % (THz) IMPORTANT: sweep in input wavelength
        % laserfreq = 192.41   % (THz) 192.41 has nice double self-pulsing
    case "newChip_SCISSOR-2-1"
        laserfreq_sweep = 192.395:0.005:192.435;   % (THz) IMPORTANT: sweep in input wavelength
    case "newChip_RING-1-3"
        laserfreq_sweep = 192.98:0.005:193.02;   % (THz) IMPORTANT: sweep in input wavelength
    case "newChip_SCISSOR-3-multiouts"
        laserfreq_sweep = 192.43:0.005:192.47;   % (THz) IMPORTANT: sweep in input wavelength
    case "newChip_RING-lowQ"
        laserfreq_sweep = 192.36:0.0025:192.44;   % (THz) IMPORTANT: sweep in input wavelength
    case "structure_1"
        laserfreq_sweep = 192.36:0.005:192.44;   % (THz) IMPORTANT: sweep in input wavelength
end   % end switch structure

VOAin_pwr_sweep = 0.1:1/15:1;   % IMPORTANT: sweep of input power
VOAout_pwr_ATTENTION = 0.05;   % 0.02 is fine once the fast amplified detector has a 50 ohm impendance on the picoscope, and the VOAin power spans between 0.1 and 0.6
pico.SetChannel(enable = [1,1,1,1], range = ["5V", "5V", "5V", "5V"]);
pico.SetChannel(enable = [1,0,1,1]);
laser_power = 10*100;  % dBm*100

sampling_step_approx = 1e-9   % (s) IMPORTANT. 1 ns --> 800 ps
n_acquired_points = 1e5   % IMPORTANT

% AWG_data_name = "AWGdata_SPmap_suggestedSampleRate50ns_inverted"

% AWG.set_output("C1","ON");

%% Save folder
saveFolder_date = saveFolder +structure+"_" +...
    string(datetime("now","format", 'yyyy-MM-dd_HH-mm')) + "\";
[~,~] = mkdir(saveFolder_date);

%% set instrumentation
VOAin(0);   % these are set here to 0 to enable zero-input PD reading
VOAout(0);
pause(0.3);

% AWG.set_signal_samplerate( "C1", "TARB", 2e7, 'HOLD' );   % 2e7 --> bits of 50 ns
% AWG.set_base_wave_parameters("C1","AMP", 14);   % here set Vpp of AWG output (14 if you want to go to 0 power when signal is 0)
% AWG.set_output("C1","ON");

% check_awg_name = 0
% while ~check_awg_name   % here make sure the correct awg segment is generated!
%     AWG.set_arbitrary_wave( "C1", "NAME", AWG_data_name);
%     pause(0.5);
%     a=AWG.get_arbitrary_wave('C1');
%     check_awg_name = strcmp(strip(string(a(3))), AWG_data_name +".bin");
% end
% disp('OK')

% TLS.enable
TLS.pwr(laser_power);
TLS.wait_pending_operation;


%% Calculate zeroPoints of detectors

pico.set_closest_sampling_time(sampling_step_approx*100);
pico.samples_per_trace(n_acquired_points/100);

TLS.disable;
TLS.wait_pending_operation;

pause(0.5);
pico.acquire(autoRange = 1);

zeroAmplifiedDrop = pico.channels.data{1};
zeroIn = pico.channels.data{3};
zeroThorlabsThrough = pico.channels.data{4};
% disp([zeroAmplifiedDrop, zeroIn, zeroThorlabsThrough]);

TLS.enable;
TLS.wait_pending_operation;
VOAin(VOAin_pwr_sweep(1));
VOAout(VOAout_pwr_ATTENTION);
pause(0.5);



%% measurement: nested loops to span wavelength, input power 
pico.set_closest_sampling_time(sampling_step_approx);
pico.samples_per_trace(n_acquired_points);

% get input example
pico.acquire(autoRange = 1);
input_signal = pico.channels.data{3};  
% save parameters
save(saveFolder_date + "Parameters_SPmap_2outs",...
                            "gain_PDin", "zeroAmplifiedDrop", "zeroIn", "zeroThorlabsThrough", "VOAin_pwr_sweep",...
                            "laserfreq_sweep", "VOAout_pwr_ATTENTION", "sampling_step_approx", "n_acquired_points", "input_signal", "EDFA_current" );

for i_freq=1:length(laserfreq_sweep)   % laser frequency loop
    % set laser frequency
    TLS.cleanMode(0);   % 0 --> dither noisy mode: needs to be used sometimes so that the laser relocks
    TLS.wait_pending_operation;
    TLS.freq(laserfreq_sweep(i_freq));
    TLS.wait_pending_operation;
    TLS.cleanMode(2);   % 2 --> low-noise whispering cleanmode
    TLS.wait_pending_operation;

    fprintf('****** Laser frequency: %3.3f ******* \n',laserfreq_sweep(i_freq));

    for i_pow=1:length(VOAin_pwr_sweep)   % input power loop
        VOAin(0);   % pause input so as to let rings discharge
        pause(0.3)
        VOAin(VOAin_pwr_sweep(i_pow));   % set input power
        pause(0.3)
    
        pico.acquire(autoRange=true);
        out_drop = pico.channels.data{1}; 
        % ensure the sensitive PD does not saturate:
        while (max(out_drop) > saturation_voltage) || (max(out_drop) < saturation_voltage/3) || VOAout_pwr_ATTENTION == 1
            if max(out_drop) > saturation_voltage
                VOAout_pwr_ATTENTION = VOAout_pwr_ATTENTION*0.5
            elseif max(out_drop) < saturation_voltage/2
                VOAout_pwr_ATTENTION = min(VOAout_pwr_ATTENTION*1.25, 1);
            end
            VOAout(VOAout_pwr_ATTENTION);
            disp("Adjusting VOAout to fit amplified PD range: " + string(max(out_drop)))
            pause(0.5);
            pico.acquire(autoRange=true);
            out_drop = pico.channels.data{1}; 
        end   % while max(out_drop) > saturation_voltage ...
        out_through = pico.channels.data{4};
        input_average = mean(pico.channels.data{3});
        % out_add = pico.channels.data{4};
        
        save(saveFolder_date + "SPmap_2outs_iFreq"+string(i_freq)+"_iPow"+string(i_pow),...
                        "out_drop", "out_through", "input_average", "VOAout_pwr_ATTENTION");
        
    end   % for i_pow
end   % for i_freq

VOAout(0.01);
disp("*** End of measurement ***")

%%   visualize saved data
% i_freq = 8
% i_pow = 12
% load("Saved_data/1R_SPmap_class/TEST_structure_1_2023-12-12_11-26\" + "SPmap_2outs_iFreq"+string(i_freq)+"_iPow"+string(i_pow))
% hold on
% plot(out_through)
% plot(out_drop)




