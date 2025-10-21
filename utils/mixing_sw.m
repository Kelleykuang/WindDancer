function array_mix_sw = mixing_sw(array_rece_sw, trans_sw_cos, trans_sw_sin, ...
                                   paras, Vs, min_range, max_range, num_of_tabs)
% MIXING_SW FMCW demodulation via I/Q mixing
%
% DESCRIPTION:
%   Performs FMCW demodulation by mixing received acoustic signal with
%   reference chirp signals (cosine and sine). This extracts beat frequencies
%   that are proportional to target distance. The process includes bandpass
%   filtering to isolate the chirp frequency band and lowpass filtering to
%   extract beat frequencies corresponding to the valid range window.
%
% SYNTAX:
%   array_mix_sw = mixing_sw(array_rece_sw, trans_sw_cos, trans_sw_sin, ...
%                            paras, Vs, min_range, max_range, num_of_tabs)
%
% INPUTS:
%   array_rece_sw  - Received audio signal [N x 1]
%   trans_sw_cos   - Reference chirp cosine component [M x 1]
%   trans_sw_sin   - Reference chirp sine component [M x 1]
%   paras          - Structure with chirp parameters (.B, .T, .Fs, .Fc)
%   Vs             - Speed of sound [m/s] (e.g., 343)
%   min_range      - Minimum detection range [m]
%   max_range      - Maximum detection range [m]
%   num_of_tabs    - FIR filter order (e.g., 100)
%
% OUTPUTS:
%   array_mix_sw   - Complex IF (intermediate frequency) signal [M x K]
%                    M = samples per chirp, K = number of chirps
%
% PROCESSING STEPS:
%   1. Bandpass filter: Isolate chirp frequency band [Fc to Fc+B]
%   2. I/Q mixing: Multiply received signal by cos and sin references
%   3. Lowpass filter: Extract beat frequencies for valid range window
%   4. Combine I/Q: Create complex baseband signal
%
% BEAT FREQUENCY RELATIONSHIP:
%   f_beat = 2*B*R/(Vs*T)
%   where R is target range [m]
%
% EXAMPLE:
%   [cos_ref, sin_ref, ~] = generate_transmit_sw(paras);
%   if_signal = mixing_sw(audio, cos_ref, sin_ref, paras, 343, 0.1, 3, 100);
%   range_fft = fft(if_signal, [], 1);  % Range FFT
%
% SEE ALSO:
%   generate_transmit_sw, fir1, filtfilt
%
% REFERENCE:
%   - FMCW radar demodulation and beat frequency extraction
%   - I/Q (In-phase/Quadrature) signal processing

    % Extract parameters
    B = paras.B;          % Bandwidth [Hz]
    T = paras.T;          % Chirp duration [s]
    Fs = paras.Fs;        % Sampling frequency [Hz]
    Fc = paras.Fc;        % Carrier frequency [Hz]
    
    % Calculate samples per chirp
    single_chirp_len = T * Fs;
    
    %% Design Filters
    
    % Bandpass filter: Isolate chirp frequency band [Fc-10 to Fc+B+10]
    % This removes out-of-band noise and interference
    bpFilt_chirp = fir1(num_of_tabs, [(Fc-10)/(Fs/2), (Fc+B+10)/(Fs/2)], 'bandpass');
    
    % Beat frequency range corresponding to distance range
    % f_beat = 2*B*R/(Vs*T)
    low_freq = min_range * 2/Vs/T * B;   % Beat freq for min_range
    high_freq = max_range * 2/Vs/T * B;  % Beat freq for max_range
    
    % Lowpass filter: Extract beat frequencies for valid range window
    lpFilt_chirp = fir1(num_of_tabs, [(low_freq)/(Fs/2), (high_freq)/(Fs/2)], 'bandpass');
    
    %% Process Each Chirp
    
    % Calculate number of complete chirps in signal
    num_of_chirps = floor(size(array_rece_sw, 1) / single_chirp_len);
    
    % Initialize output array (complex IF signal)
    array_mix_sw = zeros(single_chirp_len, num_of_chirps);
    
    % Process each chirp
    for chirp_idx = 1:num_of_chirps-1
        % Extract single chirp from received signal
        indices = (chirp_idx-1)*single_chirp_len+1 : chirp_idx*single_chirp_len;
        chirp_array_rece_sw = array_rece_sw(indices, :);
        
        % Step 1: Bandpass filter received signal
        chirp_array_rece_sw = filtfilt(bpFilt_chirp, 1, chirp_array_rece_sw);
        
        % Step 2: I/Q mixing (multiply by reference signals)
        mix_sw_cos = chirp_array_rece_sw .* trans_sw_cos;  % In-phase component
        mix_sw_sin = chirp_array_rece_sw .* trans_sw_sin;  % Quadrature component
        
        % Step 3: Lowpass filter to extract beat frequencies
        mix_sw_cos = filtfilt(lpFilt_chirp, 1, mix_sw_cos);
        mix_sw_sin = filtfilt(lpFilt_chirp, 1, mix_sw_sin);
        
        % Step 4: Combine I/Q into complex signal
        array_mix_sw(:, chirp_idx) = mix_sw_cos + 1j*mix_sw_sin;
    end
end
