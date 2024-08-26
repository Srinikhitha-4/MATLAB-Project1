% Hearing Aid Simulation

% Parameters
Fs = 44200; % Sampling frequency
duration = 2; % Duration in seconds
delay = 0.2; % Delay for feedback (seconds)
feedback_gain = 0.5; % Feedback gain
drc_threshold_dB = -30; % DRC threshold in dB
drc_ratio = 5; % DRC ratio
noise_reduction = true; % Toggle noise reduction

% Generate a sample input signal (e.g., a sinusoidal wave)
t = 0:1/Fs:duration-1/Fs; % Time vector
input_signal = sin(2*pi*440*t); % Example input signal (440 Hz tone)

% Create a delay line for feedback
feedback_delay = round(delay * Fs); % Convert delay to samples
feedback_buffer = zeros(1, feedback_delay);

% Initialize output signal
output_signal = zeros(1, length(input_signal));
%% 
% Process signal with feedback loop
for i = 1:length(input_signal)
    % Feedback signal
    feedback_signal = feedback_gain * feedback_buffer(1);
    
    % Output signal calculation
    output_signal(i) = input_signal(i) + feedback_signal;
    
    % Update feedback buffer
    feedback_buffer = [output_signal(i), feedback_buffer(1:end-1)];
end
%% 
% Dynamic Range Compression (DRC)
% Convert DRC threshold to linear scale
drc_threshold = 10^(drc_threshold_dB / 30);

% Initialize compressed signal
compressed_signal = zeros(size(output_signal));

% Apply DRC
for i = 1:length(output_signal)
    % Check if the signal exceeds the threshold
    if abs(output_signal(i)) > drc_threshold
        % Apply compression
        compressed_signal(i) = sign(output_signal(i)) * ...
            (drc_threshold + (abs(output_signal(i)) - drc_threshold) / drc_ratio);
    else
        % No compression applied
        compressed_signal(i) = output_signal(i);
    end
end
%% 
% Noise Reduction (NR)
if noise_reduction
    % Parameters for Spectral Subtraction
    noise_estimate = 0.02; % Example noise level (adjust as needed)
    
    % Compute power spectrum of the noisy signal
    N = length(compressed_signal);
    X = fft(compressed_signal);
    P = abs(X).^2 / N; % Power spectrum
    
    % Spectral subtraction
    P_clean = max(P - noise_estimate, 0); % Subtract noise estimate and avoid negative values
    
    % Compute the inverse FFT to get the time-domain signal
    X_clean = sqrt(P_clean) .* exp(1i * angle(X));
    clean_signal = real(ifft(X_clean, 'symmetric'));
    
    % Ensure the clean signal has the same length as the input signal
    clean_signal = clean_signal(1:N);
else
    clean_signal = compressed_signal;
end
%% 
% Feedback Suppression (DFS/DFC)
% Design a notch filter for feedback suppression
f0 = 1000; % Frequency of feedback (example)
Q = 30; % Quality factor
notch_filter = designfilt('bandstopiir', 'FilterOrder', 2, 'HalfPowerFrequency1', f0-50, 'HalfPowerFrequency2', f0+50, 'DesignMethod', 'butter', 'SampleRate', Fs);
filtered_signal = filtfilt(notch_filter, clean_signal);
%% 
% Plotting the results
figure;
subplot(4,1,1);
plot(t, input_signal);
title('Input Signal');
xlabel('Time (sec)');
ylabel('Amplitude');

subplot(4,1,2);
plot(t, output_signal);
title('Output Signal with Feedback Loop');
xlabel('Time (sec)');
ylabel('Amplitude');

subplot(4,1,3);
plot(t, compressed_signal);
title('Compressed Output Signal (DRC)');
xlabel('Time (sec)');
ylabel('Amplitude');

subplot(4,1,4);
plot(t, filtered_signal);
title('Final Output Signal with Noise Reduction and Feedback Suppression');
xlabel('Time (sec)');
ylabel('Amplitude');

% Save the final output signal to a file
audiowrite('final_output_signal.wav', filtered_signal, Fs);

disp('Simulation complete. The final output signal is saved as ''final_output_signal.wav''.');






