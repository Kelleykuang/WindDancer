%% ========================================================================
% PHASE TRACKING FOR MICRO-MOTION DETECTION (VMD-BASED)
% =========================================================================
%
% DESCRIPTION:
%   Evaluates phase-based micro-motion tracking on real-world acoustic data
%   using Variational Mode Decomposition (VMD). This script demonstrates
%   advanced signal decomposition for isolating breathing or heartbeat
%   signals from complex phase measurements with multiple interference
%   components.
%
% PURPOSE:
%   - Extract micro-motion from phase-unwrapped signals
%   - Apply VMD to separate signal components (breathing vs. noise)
%   - Compare raw phase tracking vs. VMD-enhanced tracking
%   - Evaluate peak-to-peak amplitude accuracy
%
% METHODOLOGY:
%   1. Load and synchronize real-world audio data
%   2. Apply FMCW demodulation to extract IF signal
%   3. Perform range FFT and extract target bin phase
%   4. Baseline: Direct phase unwrapping
%   5. Advanced: VMD decomposition to isolate motion component
%   6. Evaluate both methods using peak detection
%
% KEY INNOVATION:
%   VMD (Variational Mode Decomposition) separates complex phase signals
%   into intrinsic mode functions (IMFs), isolating breathing motion from
%   environmental noise, phase drift, and other artifacts.
%
% INPUTS:
%   - Real-world audio recordings in 'realworld_data/'
%   - Known target location and motion frequency
%
% OUTPUTS:
%   - errors_baselines: Errors for raw phase tracking
%   - errors_methods: Errors for VMD-enhanced tracking
%   - Visualizations comparing both methods
%
% SEE ALSO:
%   VMD_complex, mixing_sw, generate_transmit_sw, fine_grained_eval
%
% REFERENCE:
%   - K. Dragomiretskiy and D. Zosso, "Variational Mode Decomposition,"
%     IEEE Trans. Signal Process., 2014.
%
% =========================================================================

%% Initialization
clear; close all; clc;
addpath('utils');  % Add utility functions to path

%% =====================================================================
%  SYSTEM PARAMETERS
%  =====================================================================

% FMCW chirp parameters (longer chirp for better phase resolution)
paras.B = 4e3;              % Bandwidth [Hz]
paras.T = 0.15;             % Chirp duration [s] (longer than ranging)
paras.Fs = 48e3;            % Sampling frequency [Hz]
paras.Fc = 18e3;            % Carrier frequency [Hz]

% Speed of sound
Vs = 343;                   % [m/s] at 20°C

% Derived parameters
chirp_len = paras.T * paras.Fs;  % Samples per chirp

%% =====================================================================
%  TARGET CONFIGURATION
%  =====================================================================

% Known target location (ground truth)
% target_loc = 0.985;       % Alternative target location [m]
target_loc = 1.37;          % Target distance [m]

% Known motion frequency (e.g., from breathing rate measurement)
freq_est = 0.729;           % Motion frequency [Hz] (~44 breaths/min)

%% =====================================================================
%  PROCESSING PARAMETERS
%  =====================================================================

% Filtering parameters
filter_tabs = 100;          % FIR filter order

% Range parameters
min_range = 0.1;            % Minimum detection range [m]
max_range = min([3, 0.5*paras.T*Vs]);  % Maximum range [m]

% Calculate distance bins
dists = (0:chirp_len-1) * Vs / (2*paras.B);
ind_range = ((dists >= min_range) & (dists <= max_range));
dists_range = dists(ind_range);

% Find target bin for phase extraction
[~, target_bin] = min(abs(dists_range - target_loc));

%% =====================================================================
%  DATA SELECTION
%  =====================================================================

% Time window for processing
start_t = 1*paras.Fs + 1;   % Skip first second
end_t = start_t + 121*paras.Fs;  % Process 121 seconds

%% =====================================================================
%  GENERATE REFERENCE SIGNALS
%  =====================================================================

% Generate transmitted chirp signals (cosine and sine for I/Q mixing)
[trans_sw_cos, trans_sw_sin, ~] = generate_transmit_sw(paras);

%% =====================================================================
%  VMD AND EVALUATION PARAMETERS
%  =====================================================================

% Parameters for signal decomposition and evaluation
segment_size = 0.5 / paras.T;    % Segment size for processing
half_move = 0.005;                % Half of peak-to-peak amplitude [m] (5mm)
numImfs = 4;                      % Number of intrinsic mode functions for VMD

%% =====================================================================
%  DATA LOADING AND PROCESSING
%  =====================================================================

% Data path and file selection
datapath = 'realworld_data/';
filenames = ["example_2.wav"];   % Real-world recording with breathing motion

% Initialize result arrays
errors_baselines = [];    % Errors for raw phase tracking
std_baselines = [];       % Standard deviations for baseline
errors_methods = [];      % Errors for VMD-enhanced method
std_methods = [];         % Standard deviations for VMD method

fprintf('Processing %d file(s) for phase-based micro-motion tracking...\n', length(filenames));

% Process each file
for i = 1:length(filenames)
    fprintf('\nProcessing file %d/%d: %s\n', i, length(filenames), filenames(i));
    
    %% -----------------------------------------------------------------
    %  LOAD AND EXTRACT AUDIO SEGMENT
    %  -----------------------------------------------------------------
    
    % Read audio file
    [array_rece_sw, ~] = audioread(strcat(datapath, filenames(i)));
    
    % Extract specified time window
    array_rece_sw = array_rece_sw(start_t:end, 1);
    
    %% -----------------------------------------------------------------
    %  SYNCHRONIZATION: FIND START OF CHIRP SEQUENCE
    %  -----------------------------------------------------------------
    
    % Use cross-correlation to find where chirp sequence starts
    [c, lags] = xcorr(array_rece_sw(1:10*chirp_len), trans_sw_cos);
    
    % Keep only positive lags
    c = c(ceil(length(lags)/2):end);
    lags = lags(ceil(length(lags)/2):end);
    
    % Find peak correlation (synchronization point)
    [~, max_idx] = max(c);
    num_of_start_delay_samples = lags(max_idx) + 1;
    
    % Align received signal with transmitted chirp
    array_rece_sw = array_rece_sw(num_of_start_delay_samples:end, :);
    
    %% -----------------------------------------------------------------
    %  FMCW DEMODULATION AND RANGE FFT
    %  -----------------------------------------------------------------
    
    % Apply I/Q mixing to extract beat frequencies
    if_signal = mixing_sw(array_rece_sw, trans_sw_cos, trans_sw_sin, ...
                          paras, Vs, min_range, max_range, filter_tabs);
    
    % Transform to range domain
    range_fft = fft(if_signal, chirp_len, 1);
    range_fft = range_fft(ind_range, :);
    
    %% -----------------------------------------------------------------
    %  BASELINE METHOD: RAW PHASE UNWRAPPING
    %  -----------------------------------------------------------------
    
    % Extract complex signal at target bin (time series)
    cfrs = range_fft(target_bin, :);
    
    % Detrend I and Q components separately
    cfrs_detrend = detrend(real(cfrs), 0) + 1j*detrend(imag(cfrs), 0);
    
    % Unwrap phase and convert to displacement
    % Displacement = λ/(4π) * Δφ, where λ = c/f
    phase = unwrap(angle(cfrs_detrend));
    movement_raw = phase * Vs / (4*pi*paras.Fc);
    
    fprintf('  Baseline: Raw phase unwrapping complete\n');
    
    %% -----------------------------------------------------------------
    %  ADVANCED METHOD: VMD-BASED SIGNAL DECOMPOSITION
    %  -----------------------------------------------------------------
    
    fprintf('  Applying VMD decomposition...\n');
    tic;
    
    % Apply Variational Mode Decomposition to complex signal
    [imfs, res] = VMD_complex(cfrs, numImfs, Vs/(4*pi*paras.Fc), ...
                              1/paras.T, 0.15);
    
    % Reconstruct signal from selected IMFs (exclude high-frequency noise)
    recon = sum(imfs(:, 1:end-1), 2);
    
    vmd_time = toc;
    fprintf('  VMD processing time: %.2f seconds\n', vmd_time);
    
    % Unwrap phase of reconstructed signal and convert to displacement
    phase = unwrap(angle(recon));
    movement = phase' * Vs / (4*pi*paras.Fc);
    
    %% -----------------------------------------------------------------
    %  EVALUATION: PEAK DETECTION AND ERROR CALCULATION
    %  -----------------------------------------------------------------
    
    % Evaluate VMD-enhanced method
    [errors] = fine_grained_eval(movement, half_move, freq_est, 1/paras.T);
    
    % Evaluate baseline method
    [errors_raw] = fine_grained_eval(movement_raw, half_move, freq_est, 1/paras.T);
    
    % Store results
    errors_baselines = [errors_baselines errors_raw];
    errors_methods = [errors_methods errors];
    
    fprintf('  Baseline mean error: %.4f mm\n', mean(errors_raw)*1000);
    fprintf('  VMD mean error: %.4f mm\n', mean(errors)*1000);
    
    %% -----------------------------------------------------------------
    %  VISUALIZATION: COMPARE METHODS
    %  -----------------------------------------------------------------
    
    figure('Position', [100, 100, 1000, 600]);
    tiledlayout(2, 1, 'TileSpacing', 'compact');
    
    % Plot 1: Movement signals comparison
    nexttile
    plot(movement_raw, 'b-', 'LineWidth', 1.5); 
    hold on;
    plot(movement, 'r-', 'LineWidth', 1.5);
    xlabel('Frame Index');
    ylabel('Displacement [m]');
    title('Micro-Motion Tracking: Raw vs. VMD-Enhanced');
    legend('Raw Phase', 'VMD-Enhanced', 'Location', 'best');
    grid on;
    
    % Plot 2: Error CDFs
    nexttile
    cdfplot(errors_raw * 1000);  % Convert to mm
    hold on;
    cdfplot(errors * 1000);
    xlabel('Error [mm]');
    ylabel('Cumulative Probability');
    title('Cumulative Distribution of Peak-to-Peak Errors');
    legend('Raw Phase', 'VMD-Enhanced', 'Location', 'best');
    grid on;
end

%% =====================================================================
%  OVERALL COMPARISON: ALL FILES
%  =====================================================================

fprintf('\n========================================\n');
fprintf('OVERALL RESULTS (NORMALIZED)\n');
fprintf('========================================\n');

% Create comparison figure
figure('Position', [100, 100, 800, 600]);
cdfplot(errors_baselines / (2*half_move)); 
hold on;
cdfplot(errors_methods / (2*half_move));
xlabel('Normalized Error (as fraction of true amplitude)');
ylabel('Cumulative Probability');
title('Overall Performance: Raw Phase vs. VMD-Enhanced');
legend('Raw Phase', 'VMD-Enhanced', 'Location', 'best');
grid on;

% Display statistics
fprintf('Baseline (Raw Phase):\n');
fprintf('  Mean Normalized Error: %.4f\n', mean(errors_baselines / (2*half_move)));
fprintf('  Mean Absolute Error: %.4f mm\n', mean(errors_baselines)*1000);
fprintf('\nVMD-Enhanced Method:\n');
fprintf('  Mean Normalized Error: %.4f\n', mean(errors_methods / (2*half_move)));
fprintf('  Mean Absolute Error: %.4f mm\n', mean(errors_methods)*1000);
fprintf('\nImprovement:\n');
fprintf('  Error Reduction: %.2f%%\n', ...
        (mean(errors_baselines) - mean(errors_methods)) / mean(errors_baselines) * 100);
fprintf('========================================\n');


%% =====================================================================
%  NOTE: HELPER FUNCTIONS
%  =====================================================================
% The following helper functions are used in this script:
%   - generate_transmit_sw.m : Generate transmitted chirp signals
%   - mixing_sw.m            : FMCW demodulation (I/Q mixing)
%   - VMD_complex.m          : Variational Mode Decomposition for complex signals
%   - fine_grained_eval.m    : Fine-grained evaluation of micro-motion detection accuracy
%
% These functions are located in the utils/ directory.
% See their individual files for detailed documentation.