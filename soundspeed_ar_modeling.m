%% ========================================================================
% SPEED-OF-SOUND AR MODEL ESTIMATION
% =========================================================================
%
% DESCRIPTION:
%   Estimates autoregressive (AR) model parameters from speed-of-sound
%   measurements collected at different conditions (distances or chirp
%   lengths). The AR models capture temporal correlations in speed-of-sound
%   variations caused by environmental factors.
%
% PURPOSE:
%   - Fit AR models to real-world speed-of-sound measurements
%   - Generate AR coefficients and noise statistics for simulation
%   - Create parameterized models for different experimental conditions
%
% METHODOLOGY:
%   1. Load raw speed-of-sound measurements
%   2. For each measurement set: detrend and fit AR model using Yule-Walker
%   3. Extract AR coefficients and noise standard deviation
%   4. Save parameters for use in simulation scripts
%
% INPUTS:
%   - Raw speed-of-sound measurements (MAT file)
%   - Expected variable: sos_est_all (cell array of measurements)
%
% OUTPUTS:
%   - ar_paras: Cell array of {AR_coefficients, noise_std}
%   - dlist or clist: List of experimental conditions (distances or chirp lengths)
%   - plist: AR model orders for each condition
%
% USAGE:
%   1. Uncomment desired dataset (distances or chirp lengths)
%   2. Run script to estimate AR parameters
%   3. Uncomment save command to save results
%
% SEE ALSO:
%   aryule, detrend, generate_ar_noise
%
% =========================================================================

%% Load Raw Measurements

% Option 1: Distance-based measurements (airflow at different distances)
load("soundspeed/raw_measurements/distances-50-230.mat");
% Expected variable: sos_est_all (cell array of speed-of-sound estimates)
dlist = (0.5:0.2:2.3);  % Distance list [m]: 0.5, 0.7, 0.9, ..., 2.3 m
plist = [5 9 9 9 1 1 1 1 1 1];  % AR model order for each distance

% Option 2: Chirp length-based measurements (different chirp durations)
% load("soundspeed/raw_measurements/chirp-lengths-30-200.mat");
% clist = [0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.20];  % Chirp durations [s]
% plist = [1 1 1 1 1 1 1 0 0];  % AR model order for each chirp length

%% AR Model Parameter Estimation

% Initialize cell array to store AR parameters
ar_paras = cell(1, length(plist));

fprintf('Estimating AR models for %d conditions...\n', length(plist));

for i = 1:length(plist)
    % Extract speed-of-sound estimates for this condition
    sos_est = sos_est_all{i};
    
    % AR model order for this condition
    p = plist(i);
    
    % Detrend signal (remove DC offset)
    signal = detrend(sos_est, 0)';
    
    % Fit AR model using Yule-Walker method
    % Returns: AR coefficients, noise variance, and reflection coefficients
    [ar_coeffs, noise_variance, ~] = aryule(signal, p);
    
    % Store AR coefficients and noise standard deviation
    ar_paras{i} = {ar_coeffs, sqrt(noise_variance)};
    
    % Display progress
    if mod(i, 5) == 0 || i == length(plist)
        fprintf('  Processed %d/%d conditions (AR order = %d)\n', i, length(plist), p);
    end
end

fprintf('AR model estimation complete!\n');

%% Display AR Model Statistics

fprintf('\nAR Model Summary:\n');
fprintf('%-10s %-12s %-15s %-15s\n', 'Condition', 'AR Order', 'Noise Std', 'First Coeff');
fprintf('%s\n', repmat('-', 1, 60));
for i = 1:length(ar_paras)
    coeffs = ar_paras{i}{1};
    noise_std = ar_paras{i}{2};
    first_coeff = coeffs(2);  % First AR coefficient (after the leading 1)
    fprintf('%-10d %-12d %-15.6f %-15.6f\n', i, plist(i), noise_std, first_coeff);
end

%% Save AR Model Parameters

% Uncomment to save results for distance-based measurements
% save("soundspeed/ar_models/distances-50-230_ar_models.mat", "ar_paras", "dlist", "plist");
% fprintf('\nSaved: soundspeed/ar_models/distances-50-230_ar_models.mat\n');

% Uncomment to save results for chirp length-based measurements
% save("soundspeed/ar_models/chirp-lengths-30-200_ar_models.mat", "ar_paras", "clist", "plist");
% fprintf('\nSaved: soundspeed/ar_models/chirp-lengths-30-200_ar_models.mat\n');

fprintf('\nNote: Uncomment save command to write results to file.\n');