function [trans_sw_cos, trans_sw_sin, trans_freq] = generate_transmit_sw(paras)
% GENERATE_TRANSMIT_SW Generate FMCW chirp reference signals
%
% DESCRIPTION:
%   Generates cosine and sine components of a linear frequency-modulated
%   chirp signal for use in FMCW (Frequency Modulated Continuous Wave)
%   acoustic ranging. The chirp sweeps linearly from carrier frequency
%   (Fc) to (Fc + B) over duration T.
%
% SYNTAX:
%   [trans_sw_cos, trans_sw_sin, trans_freq] = generate_transmit_sw(paras)
%
% INPUTS:
%   paras       - Structure containing chirp parameters:
%     .Fc       - Carrier frequency [Hz] (e.g., 18000 Hz)
%     .B        - Bandwidth [Hz] (e.g., 4000 Hz)
%     .T        - Chirp duration [s] (e.g., 0.05 s)
%     .Fs       - Sampling frequency [Hz] (e.g., 48000 Hz)
%
% OUTPUTS:
%   trans_sw_cos - Cosine component of chirp [N x 1]
%   trans_sw_sin - Sine component of chirp [N x 1]
%   trans_freq   - Instantaneous frequency over time [1 x N]
%
% CHIRP EQUATION:
%   Instantaneous frequency: f(t) = Fc + (B/T)*t
%   Phase: φ(t) = 2π*(Fc*t + B*t²/(2T))
%   Signal: s(t) = A*cos(φ(t)) + j*A*sin(φ(t))
%
% EXAMPLE:
%   paras.Fc = 18e3; paras.B = 4e3; paras.T = 0.05; paras.Fs = 48e3;
%   [cos_sig, sin_sig, freq] = generate_transmit_sw(paras);
%   plot((0:length(cos_sig)-1)/paras.Fs, cos_sig);
%
% SEE ALSO:
%   mixing_sw
%
% REFERENCE:
%   - FMCW radar signal processing
%   - Linear chirp generation

    % Signal amplitude (normalized to 1)
    amp = 1;
    
    % Initial phase offset (set to 0 for simplicity)
    init_phase = 0;
    
    % Generate time vector
    time = 0:1/paras.Fs:paras.T-1/paras.Fs;
    
    % Calculate instantaneous frequency (linear FM)
    % f(t) = Fc + (B/T)*t
    trans_freq = paras.Fc + paras.B/paras.T * time;
    
    % Calculate phase (integral of instantaneous frequency)
    % φ(t) = 2π*(Fc*t + B*t²/(2T))
    phase = paras.Fc*time + paras.B*power(time, 2)/(2*paras.T) + init_phase;
    
    % Generate cosine and sine components (I/Q signals)
    trans_sw_cos = amp * cos(2*pi*phase).';
    trans_sw_sin = amp * sin(2*pi*phase).';
end
