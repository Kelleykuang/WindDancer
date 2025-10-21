function [errors] = fine_grained_eval(movements, moving_range, f_target, fs)
% FINE_GRAINED_EVAL Evaluate micro-motion detection accuracy using peak analysis
%
% DESCRIPTION:
%   Evaluates the accuracy of phase-based micro-motion detection by
%   analyzing peaks and valleys in the reconstructed movement signal.
%   This is used for fine-grained assessment of breathing or vibration
%   detection accuracy.
%
% SYNTAX:
%   errors = fine_grained_eval(movements, moving_range, f_target, fs)
%
% INPUTS:
%   movements     - Reconstructed movement signal from phase unwrapping [1 x N]
%   moving_range  - Expected peak amplitude of movement [m]
%   f_target      - Expected frequency of movement [Hz]
%   fs            - Sampling frequency (chirp rate) [Hz]
%
% OUTPUTS:
%   errors        - Peak-to-peak amplitude errors for each cycle [1 x M]
%
% METHOD:
%   1. Detrend movement signal to remove DC offset
%   2. Detect peaks (maxima) in movement signal
%   3. Detect valleys (minima) by finding peaks of inverted signal
%   4. Calculate peak-to-peak amplitude for each cycle
%   5. Compare with expected amplitude (2*moving_range)
%
% VISUALIZATION:
%   Creates a figure with two subplots showing detected peaks and valleys
%
% EXAMPLE:
%   % Simulate breathing motion at 0.3 Hz
%   t = (0:999) / 20;  % 20 Hz sampling
%   movements = 0.005 * sin(2*pi*0.3*t);  % 5mm amplitude
%   errors = fine_grained_eval(movements, 0.005, 0.3, 20);
%   fprintf('Mean error: %.4f mm\n', mean(errors)*1000);
%
% NOTE:
%   This function is currently optional and not used by default in the
%   main simulation. Uncomment the relevant section to enable phase-based
%   evaluation.
%
% SEE ALSO:
%   generate_ar_noise, findpeaks, detrend
%
% REFERENCES:
%   - Phase-based motion estimation in FMCW radar systems
%   - Contactless vital sign monitoring using acoustic sensing

    signal_len = length(movements);
    
    % Remove DC offset
    movements = detrend(movements, 0);
    
    % Generate ground truth for reference (not currently used)
    ground_truth = moving_range .* sin(2*pi*f_target/fs*(1:signal_len*2));
    
    % Configure peak detection parameters
    npeaks = floor(f_target * 60);            % Expected number of peaks in signal
    period = fs / f_target;             % Period in samples
    min_peak_distance = period * 0.75;  % Minimum distance between peaks
    min_prominence = 0.5 * moving_range; % Minimum peak prominence
    
    % Visualize detected peaks and valleys
    figure('Name', 'Fine-Grained Movement Evaluation');
    subplot(2, 1, 1)
    findpeaks(movements, 'MinPeakDistance', min_peak_distance, ...
              'MinPeakProminence', min_prominence, 'NPeaks', npeaks);
    title('Detected Peaks (Maxima)');
    ylabel('Movement [m]');
    grid on;
    
    subplot(2, 1, 2)
    findpeaks(-movements, 'MinPeakDistance', min_peak_distance, ...
              'MinPeakProminence', min_prominence, 'NPeaks', npeaks);
    title('Detected Valleys (Minima)');
    ylabel('Movement [m]');
    xlabel('Sample Index');
    grid on;
    
    % Detect peaks (maxima)
    [peaks, ~] = findpeaks(movements, 'MinPeakDistance', min_peak_distance, ...
                          'MinPeakProminence', min_prominence, 'NPeaks', npeaks);
    
    % Detect valleys (minima)
    [valleys, ~] = findpeaks(-movements, 'MinPeakDistance', min_peak_distance, ...
                             'MinPeakProminence', min_prominence, 'NPeaks', npeaks);
    
    % Calculate peak-to-peak amplitudes
    num_peaks = length(peaks);
    num_valleys = length(valleys);
    num_pairs = min([num_peaks num_valleys]);
    
    % Error = deviation from expected peak-to-peak amplitude
    errors = abs((peaks(1:num_pairs) + valleys(1:num_pairs)) - moving_range*2);
    
    % Display mean error
    mean_error = mean(errors);
    fprintf('Mean peak-to-peak error: %.6f m (%.3f mm)\n', mean_error, mean_error*1000);

end

