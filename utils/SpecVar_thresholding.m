function [cfrs_masks, spec_vars] = SpecVar_thresholding(cfrs, segment_size, fs, windowsize, threshold)
% SPECVAR_THRESHOLDING Spectral variance-based thresholding for static clutter removal
%
% DESCRIPTION:
%   Applies spectral variance analysis to identify and mask static reflectors
%   in range-FFT data. Static objects exhibit low spectral variance, while
%   moving targets show higher variance due to Doppler shifts and micro-motion.
%   This advanced method improves upon simple background subtraction.
%
% SYNTAX:
%   [cfrs_masks, spec_vars] = SpecVar_thresholding(cfrs, segment_size, fs, ...
%                                                    windowsize, threshold)
%
% INPUTS:
%   cfrs         - Complex range-FFT data [nbins x timesteps]
%   segment_size - Number of frames per segment for analysis
%   fs           - Frame rate (chirp rate) [Hz]
%   windowsize   - Window size for spectral estimation (currently unused)
%   threshold    - Spectral variance threshold (e.g., 0.1)
%                  Lower values indicate static reflectors
%
% OUTPUTS:
%   cfrs_masks   - Binary mask [nbins x timesteps]
%                  1 = moving target (keep), 0 = static clutter (remove)
%   spec_vars    - Spectral variance metric for each bin and segment
%                  [nbins x nsegments]
%
% ALGORITHM:
%   For each range bin and time segment:
%   1. Extract time series and remove mean
%   2. Compute power spectral density (PSD) using Welch's method
%   3. Calculate spectral variance metric: var(PSD) / mean(PSD)
%   4. Compare to threshold: low variance = static, high variance = moving
%   5. Create binary mask based on threshold
%
% EXAMPLE:
%   segment_size = 10;  % 10 frames per segment
%   fs = 20;            % 20 chirps/second
%   threshold = 0.1;    % Variance threshold
%   [masks, vars] = SpecVar_thresholding(range_fft, segment_size, fs, [], threshold);
%
% SEE ALSO:
%   pwelch, var, mean
%
% REFERENCE:
%   - Spectral variance for clutter rejection in FMCW radar
%   - Advanced background subtraction techniques

    % Get dimensions
    nbins = size(cfrs, 1);      % Number of range bins
    timesteps = size(cfrs, 2);  % Number of time steps
    
    % Segment the data
    nseg = floor(timesteps / segment_size);
    cfrs_segs = reshape(cfrs(:, 1:nseg*segment_size), [nbins, nseg, segment_size]);
    
    % Initialize outputs
    cfrs_masks = ones(nbins, nseg, segment_size);  % Default: keep all data
    spec_vars = zeros(nbins, nseg);
    
    fprintf('Processing %d segments with %d range bins...\n', nseg, nbins);
    
    % Process each segment
    for s = 1:nseg
        tic;
        
        % Process each range bin
        for i = 1:nbins
            % Extract time series for this bin and segment
            seg = squeeze(cfrs_segs(i, s, :));
            
            % Remove DC component (mean)
            seg = seg - mean(seg);
            
            % Estimate power spectral density using Welch's method
            pxx = pwelch(seg, segment_size, [], segment_size);
            
            % Calculate spectral variance metric
            mean_spec = mean(pxx);
            metrics = var(pxx) / mean_spec;
            spec_vars(i, s) = metrics;
            
            % Apply threshold: low variance = static clutter
            if metrics < threshold
                cfrs_masks(i, s, :) = 0;  % Mask out static reflector
            end
        end
        
        elapsed = toc;
        if mod(s, 10) == 0 || s == nseg
            fprintf('  Segment %d/%d processed (%.2f s)\n', s, nseg, elapsed);
        end
    end
    
    % Reshape masks back to 2D
    cfrs_masks = reshape(cfrs_masks, [nbins, nseg*segment_size]);
    
    fprintf('Spectral variance thresholding complete.\n');
end


