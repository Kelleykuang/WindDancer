function noise = generate_ar_noise(coeffs, noise_std, n)
% GENERATE_AR_NOISE Generate autoregressive noise sequence
%
% DESCRIPTION:
%   Generates a time series following an autoregressive (AR) model.
%   The AR model captures temporal correlations in the noise, which
%   represents real-world speed-of-sound variations due to environmental
%   factors (temperature fluctuations, air currents, etc.).
%
% SYNTAX:
%   noise = generate_ar_noise(coeffs, noise_std, n)
%
% INPUTS:
%   coeffs    - AR model coefficients [1, a1, a2, ..., ap]
%               First coefficient should be 1 (standard AR notation)
%   noise_std - Standard deviation of driving white noise
%   n         - Length of noise sequence to generate
%
% OUTPUTS:
%   noise     - Generated AR noise sequence [n x 1]
%
% AR MODEL:
%   x(t) = -a1*x(t-1) - a2*x(t-2) - ... - ap*x(t-p) + w(t)
%   where w(t) ~ N(0, noise_std^2)
%
% EXAMPLE:
%   coeffs = [1, -0.5];  % AR(1) model with coefficient 0.5
%   noise_std = 0.1;
%   n = 1000;
%   ar_noise = generate_ar_noise(coeffs, noise_std, n);
%
% SEE ALSO:
%   fine_grained_eval
%
% REFERENCES:
%   - Box, G. E., Jenkins, G. M., Reinsel, G. C., & Ljung, G. M. (2015).
%     Time series analysis: forecasting and control. John Wiley & Sons.

    % Initialize output
    noise = zeros(n, 1);
    p = length(coeffs) - 1;  % Order of AR model
    
    if p == 0
        % Zero-order AR (white noise)
        noise = randn(n, 1) * noise_std;
    else
        % Initialize first p samples with white noise
        noise(1:p) = randn([p, 1]) * noise_std;
        
        % Generate remaining samples using AR recursion
        for i = (p+1):n
            % AR recursion: current sample depends on past p samples
            noise(i) = -coeffs(2:end)' * noise(i-p:i-1) + randn()*noise_std;
        end
    end
end

