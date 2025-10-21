%% ========================================================================
% REAL-WORLD RANGE EVALUATION WITH SPECTRAL VARIANCE FILTERING
% =========================================================================
%
% DESCRIPTION:
%   Evaluates ranging accuracy on real-world acoustic data using spectral
%   variance-based clutter rejection. This script demonstrates the complete
%   processing pipeline from raw audio to range estimation with advanced
%   static reflector suppression.
%
% PURPOSE:
%   - Process real-world acoustic recordings
%   - Apply spectral variance filtering for clutter rejection
%   - Compare baseline (simple background subtraction) vs. advanced method
%   - Evaluate ranging accuracy
%
% METHODOLOGY:
%   1. Load and preprocess recorded audio
%   2. Synchronize with transmitted chirp signal
%   3. Apply FMCW demodulation (mixing)
%   4. Perform range FFT
%   5. Apply baseline background subtraction
%   6. Apply spectral variance-based filtering
%   7. Detect target and calculate errors
%
% INPUTS:
%   - Real-world audio recordings in 'realworld_data/'
%   - Known target location for accuracy evaluation
%
% OUTPUTS:
%   - errors_baselines: Mean errors for baseline method
%   - errors_methods: Mean errors for spectral variance method
%   - Visualization of range-time maps
%
% SEE ALSO:
%   mixing_sw, generate_transmit_sw, SpecVar_thresholding
%
% =========================================================================

%% Initialization
clear; close all; clc;
addpath('utils');  % Add utility functions to path

%% =====================================================================
%  SYSTEM PARAMETERS
%  =====================================================================

% FMCW chirp parameters
paras.B = 4e3;              % Bandwidth [Hz]
paras.T = 0.05;             % Chirp duration [s]
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
target_loc = 0.985;         % Target distance [m]
% target_loc = 1.37;        % Alternative target location

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

% Find target bin (for reference)
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
%  SPECTRAL VARIANCE PARAMETERS
%  =====================================================================

% Parameters for spectral variance-based filtering
windowsize = 0.5 / paras.T;      % Window size for spectral estimation
thr = 0.1;                        % Spectral variance threshold
segment_size = 0.5 / paras.T;    % Segment size for processing

%% =====================================================================
%  DATA LOADING AND PROCESSING
%  =====================================================================

% Data path and file selection
datapath = 'realworld_data/';
filenames = ["example_1.wav"];   % Real-world recording file

% Initialize result arrays
errors_baselines = [];    % Errors for baseline method
std_baselines = [];       % Standard deviations for baseline
errors_methods = [];      % Errors for spectral variance method
std_methods = [];         % Standard deviations for spectral variance

fprintf('Processing %d file(s)...\n', length(filenames));

% Process each file
for i = 1:length(filenames)
    fprintf('\nProcessing file %d/%d: %s\n', i, length(filenames), filenames(i));
    
    %% -----------------------------------------------------------------
    %  LOAD AND EXTRACT AUDIO SEGMENT
    %  -----------------------------------------------------------------
    
    % Read audio file
    [array_rece_sw, ~] = audioread(strcat(datapath, filenames(i)));
    
    % Extract specified time window (skip initial transient)
    array_rece_sw = array_rece_sw(start_t:end_t, 1);
    
    %% -----------------------------------------------------------------
    %  SYNCHRONIZATION: FIND START OF CHIRP SEQUENCE
    %  -----------------------------------------------------------------
    
    % Use cross-correlation to find where chirp sequence starts
    [c, lags] = xcorr(array_rece_sw(1:10*chirp_len), trans_sw_cos);
    
    % Keep only positive lags (forward in time)
    c = c(ceil(length(lags)/2):end);
    lags = lags(ceil(length(lags)/2):end);
    
    % Find peak correlation (synchronization point)
    [~, max_idx] = max(c);
    num_of_start_delay_samples = lags(max_idx) + 1;
    
    % Align received signal with transmitted chirp
    array_rece_sw = array_rece_sw(num_of_start_delay_samples:end, :);
    
    %% -----------------------------------------------------------------
    %  FMCW DEMODULATION (MIXING)
    %  -----------------------------------------------------------------
    
    % Apply I/Q mixing to extract beat frequencies
    if_signal = mixing_sw(array_rece_sw, trans_sw_cos, trans_sw_sin, ...
                          paras, Vs, min_range, max_range, filter_tabs);
    
    %% -----------------------------------------------------------------
    %  RANGE FFT
    %  -----------------------------------------------------------------
    
    % Transform to range domain
    range_fft = fft(if_signal, chirp_len, 1);
    
    %% -----------------------------------------------------------------
    %  BASELINE METHOD: SIMPLE BACKGROUND SUBTRACTION
    %  -----------------------------------------------------------------
    
    % Estimate background using moving average
    mov_mean_win = 1;
    background_if_signal = movmean(if_signal, mov_mean_win, 2);
    
    % Subtract background to isolate moving target
    if_signal_diff = if_signal(:, ceil(mov_mean_win/2)+1:end) - ...
                     background_if_signal(:, 1:end-ceil(mov_mean_win/2));
    
    % Range FFT of background-subtracted signal
    range_fft_diff = fft(if_signal_diff, chirp_len, 1);
    
    % Extract valid range bins
    range_fft = range_fft(ind_range, 1:end-1);
    range_fft_diff = range_fft_diff(ind_range, 1:end-1);

    %% -----------------------------------------------------------------
    %  BASELINE: TARGET DETECTION AND ERROR CALCULATION
    %  -----------------------------------------------------------------
    
    % Calculate magnitude (dB scale for reference, not used in detection)
    magnitude = 20*log10(abs(range_fft));
    
    % Normalize magnitude for detection (peak detection on normalized data)
    magnitude_diff = abs(range_fft_diff) ./ max(abs(range_fft_diff), [], 1);
    
    % Detect target (find peak in each time step)
    [~, detected_bin] = max(magnitude_diff, [], 1);
    detected_ranges = dists_range(detected_bin);
    
    % Calculate ranging errors
    errors = abs(detected_ranges - target_loc);
    errors_baselines = [errors_baselines mean(errors)];
    std_baselines = [std_baselines std(errors)];
    
    fprintf('  Baseline method: Mean error = %.4f m, Std = %.4f m\n', ...
            mean(errors), std(errors));
    
    %% -----------------------------------------------------------------
    %  SPECTRAL VARIANCE METHOD: ADVANCED CLUTTER REJECTION
    %  -----------------------------------------------------------------
    
    % Apply spectral variance-based thresholding to identify static reflectors
    [range_fft_masks, spec_vars] = SpecVar_thresholding(range_fft, segment_size, ...
                                                         1/paras.T, windowsize, thr);
    
    % Apply mask to magnitude data
    minsize = min([size(range_fft_masks, 2), size(magnitude_diff, 2)]);
    magnitude_diff_masked = magnitude_diff(:, 1:minsize) .* range_fft_masks(:, 1:minsize);
    
    % Detect target on masked data
    [~, detected_bin] = max(magnitude_diff_masked, [], 1);
    detected_ranges = dists_range(detected_bin);
    
    % Calculate ranging errors for spectral variance method
    errors = abs(detected_ranges - target_loc);
    errors_methods = [errors_methods mean(errors)];
    std_methods = [std_methods std(errors)];
    
    fprintf('  SpecVar method: Mean error = %.4f m, Std = %.4f m\n', ...
            mean(errors), std(errors));
    
    %% -----------------------------------------------------------------
    %  VISUALIZATION
    %  -----------------------------------------------------------------
    
    figure('Position', [100, 100, 1000, 800]);
    tiledlayout(3, 1, 'TileSpacing', 'compact');
    
    % Plot 1: Baseline method range-time map
    nexttile
    imagesc(magnitude_diff);
    xlabel('Time Step');
    ylabel('Range Bin');
    title('Baseline Method: Background Subtraction');
    colorbar;
    axis xy;
    
    % Plot 2: Spectral variance method range-time map
    nexttile
    imagesc(magnitude_diff_masked);
    xlabel('Time Step');
    ylabel('Range Bin');
    title('Spectral Variance Method: Clutter-Rejected');
    colorbar;
    axis xy;
    
    % Plot 3: Spectral variance values over time
    nexttile
    imagesc(spec_vars);
    xlabel('Segment');
    ylabel('Range Bin');
    title('Spectral Variance Metric (Low = Static, High = Moving)');
    colorbar;
    axis xy;
end

%% =====================================================================
%  DISPLAY RESULTS
%  =====================================================================

fprintf('\n========================================\n');
fprintf('FINAL RESULTS\n');
fprintf('========================================\n');
fprintf('Baseline Method:\n');
fprintf('  Mean Error: %.4f m\n', mean(errors_baselines));
fprintf('  Std Error:  %.4f m\n', mean(std_baselines));
fprintf('\nSpectral Variance Method:\n');
fprintf('  Mean Error: %.4f m\n', mean(errors_methods));
fprintf('  Std Error:  %.4f m\n', mean(std_methods));
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
%   - SpecVar_thresholding.m : Spectral variance-based clutter rejection
%
% These functions are located in the utils/ directory.
% See their individual files for detailed documentation.
