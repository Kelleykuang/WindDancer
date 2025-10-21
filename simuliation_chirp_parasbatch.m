%% ========================================================================
% CHIRP-BASED ACOUSTIC SENSING SIMULATION WITH SPEED OF SOUND VARIATIONS
% =========================================================================
%
% DESCRIPTION:
%   This script simulates chirp-based acoustic sensing (FMCW - Frequency 
%   Modulated Continuous Wave) for detecting moving objects in the presence
%   of speed-of-sound variations. The speed of sound variations are modeled
%   using autoregressive (AR) parameters estimated from real-world data. 
%
% PURPOSE:
%   - Simulate acoustic ranging with time-varying speed of sound
%   - Evaluate detection accuracy under various AR noise models
%
% METHODOLOGY:
%   1. Load AR models representing speed-of-sound variations
%   2. Generate chirp signals for static reflectors and moving target
%   3. Apply AR-modeled speed-of-sound noise to received signals
%   4. Perform range FFT processing with background subtraction
%   5. Detect target location and evaluate accuracy
%
% INPUTS:
%   - ar_paras: Cell array containing AR coefficients and noise std
%   - Loaded from: "soundspeed\ar_models\distances-50-230_ar_models.mat"
%
% OUTPUTS:
%   - err_list: Mean detection errors for each AR parameter set
%   - std_list: Standard deviation of detection errors
%   - noise_std_list: Noise standard deviation for each parameter set

%% Initialization
clear all;
rng('default');  % Set random seed for reproducibility
addpath('utils');  % Add utility functions to path

%% =====================================================================
%  SIMULATION PARAMETERS
%  =====================================================================

% --- Acoustic Signal Parameters ---
c = 343;              % Speed of sound in air [m/s] at 20°C
Fs = 48e3;            % Sampling frequency [Hz]
Fc = 18e3;            % Carrier frequency [Hz] (ultrasonic range)
T = 0.05;             % Chirp duration [s] (50 ms)
B = 4e3;              % Chirp bandwidth [Hz] (4 kHz)

% --- Time and Sample Configuration ---
chirp_sample = Fs*T;           % Number of samples per chirp
num_chirps = 5*60/T;           % Total number of chirps (5 minutes worth)
t_all = (1:T*num_chirps*Fs)'/Fs;  % Complete time vector [s]

% --- Range Processing Parameters ---
dist_bin_num = 50;             % Number of distance bins
mov_mean_win = 1;              % Moving average window for background subtraction
dist_range = (0:dist_bin_num-1)*c / (2*B);  % Distance bins [m]
                                % Range resolution = c/(2*B) ≈ 4.3 cm

%% =====================================================================
%  SCENE CONFIGURATION: STATIC REFLECTORS
%  =====================================================================

% Static reflectors represent walls, furniture, or other stationary objects
static_rf_dists = [1.5 1.3 0.7 0.5];  % Distances of static reflectors [m]
static_rf_amps = [1 1 1 1];           % Relative amplitudes (all equal)
static_rf_phs = rand([length(static_rf_dists) 1]) * 2*pi;  % Random phases [rad]

%% =====================================================================
%  SCENE CONFIGURATION: MOVING TARGET
%  =====================================================================

% Moving target with sinusoidal motion (e.g., breathing, vibration)
f_target = 0.3;           % Target movement frequency [Hz] (e.g., breathing)
target_loc = 1;           % Target center location [m]
moving_range = 0.005;     % Movement amplitude [m] (5 mm peak-to-peak)
target_ph = 0;            % Initial phase [rad]
target_amp = 1;           % Signal amplitude

% Store moving object parameters
moving_object = {target_loc, moving_range, f_target};

% Generate target distance trajectory
target_dist = target_loc + moving_range*sin(2*pi*f_target*t_all);

%% =====================================================================
%  LOAD AUTOREGRESSIVE MODELS FOR SPEED-OF-SOUND VARIATIONS
%  =====================================================================

% Load pre-computed AR parameters representing speed-of-sound variations induced by airflow
load("soundspeed\ar_models\distances-50-230_ar_models.mat");
% Expected variables: 
%   - ar_paras: cell array of {coefficients, noise_std}
%   - dlist (optional): distance list corresponding to AR parameters [m]

%% =====================================================================
%  MAIN SIMULATION LOOP: EVALUATE EACH AR PARAMETER SET
%  =====================================================================

len_paras = length(ar_paras);
err_list = zeros(len_paras, 1);      % Mean detection error per AR model
std_list = zeros(len_paras, 1);      % Std deviation of detection error
movement_list = cell(len_paras, 1);  % Store movement estimates (optional)
noise_std_list = zeros(len_paras, 1);  % Store noise std for each model

fprintf('Starting simulation with %d AR parameter sets...\n', len_paras);

for ei = 1:len_paras
    
    % Extract AR coefficients and noise standard deviation
    coeffs = ar_paras{ei}{1}';
    noise_std = ar_paras{ei}{2};
    
    % Generate speed-of-sound noise using AR model
    % This represents temporal variations in sound speed
    c_noise = generate_ar_noise(coeffs, noise_std, num_chirps);
    
    %% -----------------------------------------------------------------
    %  SIGNAL GENERATION: Chirp Reception with SoS Variations
    %  -----------------------------------------------------------------
    
    signal_reve = zeros(chirp_sample, num_chirps);
    
    for ci = 1:num_chirps
        % Time indices for current chirp
        t_single_idx = ((ci-1)*chirp_sample+1:ci*chirp_sample)';
        t_single = (1:chirp_sample)'/Fs;  % Local time within chirp
        
        % --- Moving Target Signal ---
        % Beat frequency depends on distance and speed of sound
        % f_beat = 2*B*d / (c*T) where d is distance
        f_beat = 2*B*target_dist(t_single_idx) ./ ((c+c_noise(ci))*T);
        
        % Phase includes carrier and distance-dependent term
        phase_target = 2*pi*Fc*2*target_dist(t_single_idx)./(c+c_noise(ci));
        
        % Generate complex baseband signal
        signal_target = exp(1j*(2*pi*f_beat.*t_single + phase_target + target_ph));
        signal_reve(:,ci) = signal_reve(:,ci) + target_amp * signal_target;
        
        % --- Static Reflector Signals ---
        for i = 1:length(static_rf_dists)
            % Beat frequency for static reflector
            f_beat_rf = 2*B*static_rf_dists(i) ./ ((c+c_noise(ci))*T);
            
            % Phase for static reflector
            phase_rf = 2*pi*Fc*2*static_rf_dists(i)./(c+c_noise(ci));
            
            % Generate reflector signal
            signal_rf = static_rf_amps(i)*exp(1j*(2*pi*f_beat_rf.*t_single ...
                                                  + phase_rf + static_rf_phs(i)));
            signal_reve(:,ci) = signal_reve(:,ci) + signal_rf;
        end
        
        % Optional: Add white Gaussian noise
        % signal_reve(:,ci) = awgn(signal_reve(:,ci), 30, 'measured');
    end
    
    array_mix_sw = signal_reve;
    
    %% -----------------------------------------------------------------
    %  SIGNAL PROCESSING: Range FFT
    %  -----------------------------------------------------------------
    
    % Perform FFT along range dimension (each chirp)
    range_ffts = fft(array_mix_sw, chirp_sample, 1);
    range_ffts = range_ffts(1:dist_bin_num, :);
    
    %% -----------------------------------------------------------------
    %  BACKGROUND SUBTRACTION: Remove Static Clutter
    %  -----------------------------------------------------------------
    
    % Compute moving average as background estimate
    mix_sw_background = movmean(array_mix_sw, mov_mean_win, 2);
    
    % Subtract background to isolate moving target
    array_mix_sw_diff = array_mix_sw(:, ceil(mov_mean_win/2)+1:end) - ...
                        mix_sw_background(:, 1:end-ceil(mov_mean_win/2));
    
    % FFT of background-subtracted signal
    range_ffts_diff = fft(array_mix_sw_diff, chirp_sample, 1);
    range_ffts_diff = range_ffts_diff(1:dist_bin_num, :);
    
    %% -----------------------------------------------------------------
    %  TARGET DETECTION: Coarse-Grained Range Estimation
    %  -----------------------------------------------------------------
    
    % Find peak in each range profile (corresponds to target distance)
    [~, detected_bin] = max(abs(range_ffts_diff), [], 1);
    detected_ranges = dist_range(detected_bin);
    
    % Calculate detection errors
    errors = abs(detected_ranges - target_loc);
    err_list(ei) = mean(errors);
    std_list(ei) = std(errors);
    noise_std_list(ei) = noise_std;
    
    % Progress indicator
    if mod(ei, 10) == 0
        fprintf('Completed %d/%d parameter sets...\n', ei, len_paras);
    end
    
    %% -----------------------------------------------------------------
    %  ALTERNATIVE: Fine-Grained Phase-Based Movement Tracking
    %  -----------------------------------------------------------------
    % Uncomment below for phase-based micro-motion estimation
    % Note: Need to create a strong reflector located at a distance closed to the target to see the effects
    % 
%     [~, target_bin] = min(abs(dist_range - target_loc));
%     cfrs = detrend(real(range_ffts(target_bin,:)),0) + ...
%            1j*detrend(imag(range_ffts(target_bin,:)),0);
%     phase = unwrap(angle(cfrs));
%     movement = phase * c / (4*pi*(Fc+2e3));
%     errors = fine_grained_eval(movement, moving_range, f_target, 1/T);
%     err_list(ei) = mean(errors);
%     std_list(ei) = std(errors);
%     movement_list{ei} = movement;
end

fprintf('Simulation complete!\n');

%% =====================================================================
%  RESULTS VISUALIZATION
%  =====================================================================

% Create or use distance list
if ~exist('dlist', 'var')
    % If not loaded from file, create distance list (50-230 cm)
    dlist = linspace(0.5, 2.3, len_paras);  % Distance in meters
end

% Create figure
fig = figure('Position', [100, 100, 800, 600]);

% Main plot: Ranging Error vs Distance
yyaxis left
b = bar(dlist, err_list, 'FaceColor', [0.2 0.6 1]);
hold on;
ylabel('Ranging Error (m)');
ax1 = gca;
ax1.YColor = [0.2 0.6 1];

% Secondary axis: AR Noise Standard Deviation
yyaxis right
plot(dlist, noise_std_list, '-o', 'Color', [1 0.4 0], 'LineWidth', 2, ...
     'MarkerFaceColor', [1 0.4 0], 'MarkerSize', 6);
ylabel('\bf\sigma\rm (AR Noise Std)');
ax2 = gca;
ax2.YColor = [1 0.4 0];
ylim([0.01 0.05]);

% Formatting
xlabel('Distance (m)');
title('Ranging Error vs. Distance from Airflow Source');
grid on;
legend('Ranging Error', '\sigma (AR Noise)', 'Location', 'northeast');
legend boxoff;

% Figure properties for export
fontsize(fig, 11, "points");
width = 16; height = 10;
set(gcf, 'unit', 'centimeters', 'position', [10 10 width height]);

%% =====================================================================
%  SAVE RESULTS (Optional)
%  =====================================================================
% Uncomment to save results
% save('simulation_results.mat', 'dlist', 'noise_std_list', 'std_list', 'err_list', ...
%      'static_rf_dists', 'static_rf_amps', 'moving_object');

% Optional: Save figure
% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf, '-dpng', 'ranging_error_vs_distance.png', '-r600');

%% =====================================================================
%  NOTE: HELPER FUNCTIONS
%  =====================================================================
% The following helper functions are used in this script:
%   - generate_ar_noise.m  : Generate autoregressive noise sequence
%   - fine_grained_eval.m  : Evaluate micro-motion detection accuracy
%
% These functions are located in separate files in the same directory.
% See their individual files for detailed documentation.