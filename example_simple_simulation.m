%% ========================================================================
% SIMPLE EXAMPLE: CHIRP-BASED ACOUSTIC SENSING
% =========================================================================
%
% DESCRIPTION:
%   A simplified example demonstrating the basic principles of chirp-based
%   acoustic sensing without AR-modeled speed-of-sound variations. This is
%   ideal for learning the fundamentals before running the full simulation.
%
% USAGE:
%   Simply run this script in MATLAB. No additional data files required.
%
% =========================================================================

%% Clear workspace and set random seed
clear all; close all; clc;
rng(42);  % For reproducibility

fprintf('====================================================\n');
fprintf('Simple Chirp-Based Acoustic Sensing Demonstration\n');
fprintf('====================================================\n\n');

%% =====================================================================
%  CONFIGURATION
%  =====================================================================

% Acoustic parameters
c = 343;              % Speed of sound [m/s]
Fs = 48e3;            % Sampling frequency [Hz]
Fc = 18e3;            % Carrier frequency [Hz]
T = 0.05;             % Chirp duration [s]
B = 4e3;              % Bandwidth [Hz]

% Simulation parameters
chirp_sample = Fs * T;
num_chirps = 100;     % Shorter simulation for quick demo

% Range parameters
dist_bin_num = 50;
dist_range = (0:dist_bin_num-1) * c / (2*B);

fprintf('System Parameters:\n');
fprintf('  Sampling Frequency: %.1f kHz\n', Fs/1e3);
fprintf('  Carrier Frequency: %.1f kHz\n', Fc/1e3);
fprintf('  Bandwidth: %.1f kHz\n', B/1e3);
fprintf('  Chirp Duration: %.1f ms\n', T*1e3);
fprintf('  Range Resolution: %.2f cm\n', c/(2*B)*100);
fprintf('  Maximum Range: %.2f m\n\n', max(dist_range));

%% =====================================================================
%  SCENE SETUP
%  =====================================================================

% Static reflectors (walls, furniture)
static_dists = [1.5, 1.0, 0.6];
static_amps = [0.8, 1.0, 0.5];
static_phases = rand(size(static_dists)) * 2*pi;

% Moving target (simulating breathing)
target_loc = 1.2;         % Center location [m]
breathing_freq = 0.3;     % Breathing frequency [Hz] = 18 breaths/min
breathing_amp = 0.005;    % Amplitude [m] = 5 mm
target_amp = 1.0;

fprintf('Scene Configuration:\n');
fprintf('  Static reflectors at: [%.2f, %.2f, %.2f] m\n', static_dists);
fprintf('  Moving target at: %.2f m\n', target_loc);
fprintf('  Breathing frequency: %.2f Hz (%.0f breaths/min)\n', ...
        breathing_freq, breathing_freq*60);
fprintf('  Breathing amplitude: %.1f mm\n\n', breathing_amp*1000);

%% =====================================================================
%  SIGNAL GENERATION
%  =====================================================================

fprintf('Generating signals...\n');
signal_received = zeros(chirp_sample, num_chirps);

% Time vector for movement
t_all = (0:num_chirps-1) * T;
target_trajectory = target_loc + breathing_amp * sin(2*pi*breathing_freq*t_all);

for ci = 1:num_chirps
    t_single = (0:chirp_sample-1)' / Fs;
    target_dist = target_trajectory(ci);
    
    % Moving target signal
    f_beat_target = 2*B*target_dist / (c*T);
    phase_target = 2*pi*Fc*2*target_dist / c;
    signal_target = target_amp * exp(1j*(2*pi*f_beat_target*t_single + phase_target));
    signal_received(:,ci) = signal_target;
    
    % Static reflectors
    for i = 1:length(static_dists)
        f_beat_static = 2*B*static_dists(i) / (c*T);
        phase_static = 2*pi*Fc*2*static_dists(i) / c;
        signal_static = static_amps(i) * exp(1j*(2*pi*f_beat_static*t_single + ...
                                                   phase_static + static_phases(i)));
        signal_received(:,ci) = signal_received(:,ci) + signal_static;
    end
end

fprintf('Signal generation complete.\n\n');

%% =====================================================================
%  SIGNAL PROCESSING
%  =====================================================================

fprintf('Processing signals...\n');

% Range FFT (detects all objects)
range_fft = fft(signal_received, chirp_sample, 1);
range_fft = range_fft(1:dist_bin_num, :);

% Background subtraction (isolates moving target)
signal_background = movmean(signal_received, 5, 2);
signal_diff = signal_received(:, 3:end) - signal_background(:, 1:end-2);
range_fft_diff = fft(signal_diff, chirp_sample, 1);
range_fft_diff = range_fft_diff(1:dist_bin_num, :);

fprintf('Processing complete.\n\n');

%% =====================================================================
%  DETECTION
%  =====================================================================

fprintf('Detecting target...\n');

% Detect target location in each chirp
[~, detected_bins] = max(abs(range_fft_diff), [], 1);
detected_ranges = dist_range(detected_bins);

% Calculate detection accuracy
mean_detected = mean(detected_ranges);
error = abs(mean_detected - target_loc);

fprintf('Target Detection Results:\n');
fprintf('  True location: %.3f m\n', target_loc);
fprintf('  Detected location: %.3f m\n', mean_detected);
fprintf('  Detection error: %.3f cm\n\n', error*100);

%% =====================================================================
%  VISUALIZATION
%  =====================================================================

fprintf('Creating visualizations...\n\n');

% Figure 1: Range-Time Map (with all objects)
figure('Position', [100, 100, 1200, 800], 'Name', 'Acoustic Sensing Results');

subplot(2, 2, 1)
imagesc(t_all, dist_range, abs(range_fft));
xlabel('Time [s]');
ylabel('Range [m]');
title('Range-Time Map (All Objects)');
colorbar;
colormap('jet');
hold on;
plot(t_all, target_trajectory, 'w--', 'LineWidth', 2);
for i = 1:length(static_dists)
    yline(static_dists(i), 'w:', 'LineWidth', 1.5);
end
legend('Moving Target', 'Static Reflectors', 'Location', 'best');

% Figure 2: Range-Time Map (background subtracted)
subplot(2, 2, 2)
imagesc(t_all(1:end-2), dist_range, abs(range_fft_diff));
xlabel('Time [s]');
ylabel('Range [m]');
title('Range-Time Map (Moving Target Only)');
colorbar;
colormap('jet');
hold on;
plot(t_all, target_trajectory, 'w--', 'LineWidth', 2);
legend('Moving Target', 'Location', 'best');

% Figure 3: Single Range Profile
subplot(2, 2, 3)
single_profile = abs(range_fft(:, 50));
plot(dist_range, single_profile, 'b-', 'LineWidth', 2);
xlabel('Range [m]');
ylabel('Signal Magnitude');
title('Single Range Profile (Chirp #50)');
grid on;
hold on;
xline(target_loc, 'r--', 'LineWidth', 2);
for i = 1:length(static_dists)
    xline(static_dists(i), 'k:', 'LineWidth', 1.5);
end
legend('Range Profile', 'Moving Target', 'Static Reflectors', 'Location', 'best');

% Figure 4: Detected Range Over Time
subplot(2, 2, 4)
plot(t_all(1:length(detected_ranges)), detected_ranges, 'bo-', 'MarkerSize', 4);
hold on;
plot(t_all, target_trajectory, 'r--', 'LineWidth', 2);
yline(target_loc, 'k:', 'LineWidth', 1);
xlabel('Time [s]');
ylabel('Range [m]');
title('Target Detection Over Time');
legend('Detected', 'True Trajectory', 'True Center', 'Location', 'best');
grid on;
ylim([target_loc - 2*breathing_amp, target_loc + 2*breathing_amp]);

%% =====================================================================
%  MICRO-MOTION EXTRACTION (BONUS)
%  =====================================================================

% Extract phase at target bin for micro-motion
[~, target_bin] = min(abs(dist_range - target_loc));
target_signal = range_fft(target_bin, :);

% Unwrap phase to extract movement
phase_unwrapped = unwrap(angle(target_signal));
micro_motion = phase_unwrapped * c / (4*pi*Fc);
micro_motion = detrend(micro_motion);  % Remove trend

% Figure 5: Micro-motion extraction
figure('Position', [150, 150, 800, 600], 'Name', 'Micro-Motion Analysis');

subplot(2, 1, 1)
plot(t_all, target_trajectory - target_loc, 'r-', 'LineWidth', 2);
hold on;
plot(t_all, micro_motion, 'b--', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Displacement [m]');
title('Extracted Micro-Motion (Breathing Pattern)');
legend('True Motion', 'Extracted Motion', 'Location', 'best');
grid on;

subplot(2, 1, 2)
true_motion = breathing_amp * sin(2*pi*breathing_freq*t_all);
error_motion = abs(micro_motion - true_motion);
plot(t_all, error_motion * 1000, 'k-', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Error [mm]');
title('Micro-Motion Extraction Error');
grid on;

fprintf('Mean micro-motion error: %.3f mm\n', mean(error_motion)*1000);

fprintf('\n====================================================\n');
fprintf('Demonstration Complete!\n');
fprintf('====================================================\n');

