function [imfs, res] = VMD_complex(sig, numImfs, coeffs, fs, fres)
% VMD_COMPLEX Variational Mode Decomposition for complex signals
%
% DESCRIPTION:
%   Applies Variational Mode Decomposition (VMD) to a complex-valued signal
%   by separately decomposing the real and imaginary components, then
%   recombining them. This is useful for processing FMCW radar/sonar signals
%   where both I (in-phase) and Q (quadrature) components carry information.
%
% SYNTAX:
%   [imfs, res] = VMD_complex(sig, numImfs, coeffs, fs, fres)
%
% INPUTS:
%   sig      - Complex input signal [N x 1]
%   numImfs  - Number of intrinsic mode functions to extract
%   coeffs   - Conversion coefficient for phase-to-displacement
%              (typically Vs/(4*pi*Fc) for acoustic ranging)
%   fs       - Sampling frequency [Hz]
%   fres     - Frequency resolution for spectrum analysis (unused in current version)
%
% OUTPUTS:
%   imfs     - Intrinsic mode functions [N x numImfs] (complex)
%              Each column is a decomposed component
%   res      - Residual signal [N x 1] (complex)
%              Remaining signal after extracting all IMFs
%
% METHOD:
%   1. Separate complex signal into real and imaginary parts
%   2. Apply VMD independently to each part
%   3. Recombine decomposed components: IMF = Real_IMF + j*Imag_IMF
%
% VMD CONCEPT:
%   VMD decomposes a signal into narrow-band modes (IMFs) by solving
%   an optimization problem. Unlike EMD, VMD is non-recursive and more
%   robust to noise. Each IMF represents a distinct oscillatory mode.
%
% TYPICAL USAGE FOR MICRO-MOTION DETECTION:
%   - IMF 1-3: Low-frequency components (breathing, slow drift)
%   - IMF 4+: High-frequency noise, harmonics
%   - Residual: DC component or very low-frequency trend
%
% EXAMPLE:
%   sig = complex_radar_signal;  % Complex IF signal from target bin
%   numImfs = 4;                 % Decompose into 4 modes
%   coeffs = 343/(4*pi*18e3);    % Phase-to-displacement conversion
%   [imfs, res] = VMD_complex(sig, numImfs, coeffs, 20, 0.15);
%   breathing_signal = sum(imfs(:,1:2), 2);  % Combine first 2 IMFs
%
% ALTERNATIVE METHODS (COMMENTED OUT):
%   - EMD (Empirical Mode Decomposition): Classic adaptive decomposition
%   - CEEMDAN: Noise-assisted EMD variant
%   - EWT (Empirical Wavelet Transform): Wavelet-based decomposition
%
% SEE ALSO:
%   vmd, emd, findpeaks, unwrap
%
% REFERENCE:
%   - K. Dragomiretskiy and D. Zosso, "Variational Mode Decomposition,"
%     IEEE Trans. Signal Process., vol. 62, no. 3, pp. 531-544, 2014.
%
% NOTE:
%   Requires Signal Processing Toolbox for vmd() function (R2021a+)

    %% Apply VMD to Real and Imaginary Components Separately
    
    % Decompose real part
    [real_imfs, real_res, ~] = vmd(real(sig), "NumIMFs", numImfs);
    
    % Decompose imaginary part
    [imag_imfs, imag_res, ~] = vmd(imag(sig), "NumIMFs", numImfs);
    
    %% Alternative Decomposition Methods (Commented Out)
    
    % Option 1: Empirical Mode Decomposition (EMD)
    % [real_imfs, real_res] = emd(real(sig), "MaxNumIMF", numImfs);
    % [imag_imfs, imag_res] = emd(imag(sig), "MaxNumIMF", numImfs);
    
    % Option 2: Complete Ensemble EMD with Adaptive Noise (CEEMDAN)
    % Nstd = 0.1;      % Noise standard deviation
    % NR = 50;         % Number of realizations
    % MaxIter = 300;   % Maximum iterations
    % SNRFlag = 1;     % SNR flag
    % [real_imfs, ~] = ceemdan(real(sig), Nstd, NR, MaxIter, SNRFlag);
    % [imag_imfs, ~] = ceemdan(imag(sig), Nstd, NR, MaxIter, SNRFlag);
    
    % Option 3: Empirical Wavelet Transform (EWT)
    % [imfs, cfs] = ewt(sig, 'MaxNumPeaks', numImfs);
    
    %% Recombine Complex IMFs
    
    % Combine real and imaginary IMFs
    imfs = real_imfs + 1j*imag_imfs;
    
    % Combine residuals
    res = real_res + 1j*imag_res;
    
    %% Optional Visualization (Commented Out)
    % Uncomment below to visualize IMF decomposition and spectra
    %
    % freq_range = [-150 0];
    % figure('Name', 'VMD Decomposition Analysis');
    % tiledlayout(size(imfs, 2)+2, 2);
    % 
    % % Original signal spectrum and phase
    % nexttile
    % pspectrum(sig, fs, 'FrequencyResolution', fres);
    % ylim(freq_range);
    % title('Original Signal Spectrum');
    % 
    % nexttile
    % plot(coeffs*unwrap(angle(sig)));
    % title('Original Signal Phase');
    % ylabel('Displacement [m]');
    % 
    % % Each IMF's spectrum and phase
    % for i = 1:size(imfs, 2)
    %     nexttile
    %     pspectrum(imfs(:,i), fs, 'FrequencyResolution', fres);
    %     ylim(freq_range);
    %     title(sprintf('IMF %d Spectrum', i));
    %     
    %     nexttile
    %     plot(coeffs*unwrap(angle(imfs(:,i))));
    %     title(sprintf('IMF %d Phase', i));
    %     ylabel('Displacement [m]');
    % end
    % 
    % % Residual spectrum and phase
    % nexttile
    % pspectrum(res, fs);
    % ylim(freq_range);
    % title('Residual Spectrum');
    % 
    % nexttile
    % plot(coeffs*unwrap(angle(res)));
    % title('Residual Phase');
    % ylabel('Displacement [m]');
end
