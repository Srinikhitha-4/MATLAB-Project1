% Hearing Aid Simulation with Feedback Loop, DRC, Noise Reduction, and Feedback Suppression

% Parameters
Fs = 44200; % Sampling frequency
duration = 2; % Duration in seconds
delay = 0.01; % Delay for feedback (in seconds) reduced for realistic hearing aid
feedback_gain = 0.5; % Feedback gain
drc_threshold_dB = -30; % DRC threshold in dB
drc_ratio = 5; % DRC ratio
attack_time = 0.01; % Compressor attack time in seconds
release_time = 0.1; % Compressor release time in seconds
noise_reduction = true; % Toggle noise reduction

% Generate a sample input signal (e.g., a sinusoidal wave)
t = 0:1/Fs:duration-1/Fs; % Time vector
input_signal = sin(2*pi*440*t); % Example input signal (440 Hz tone)

% Create a delay line for feedback
feedback_delay = round(delay * Fs); % Convert delay to samples
feedback_buffer = zeros(1, feedback_delay);

% Initialize output signal
output_signal = zeros(1, length(input_signal));
compressed_signal = zeros(1, length(input_signal));

% Process signal with feedback loop
for i = 1:length(input_signal)
    % Feedback signal
    feedback_signal = feedback_gain * feedback_buffer(1);
    
    % Output signal with feedback
    output_signal(i) = input_signal(i) + feedback_signal;
    
    % Dynamic Range Compression
    signal_rms = sqrt(mean(output_signal(max(1, i-10):i).^2)); % RMS windowing
    if signal_rms > 10^(drc_threshold_dB/20)
        compression_gain = (10^(drc_threshold_dB/20)) + (signal_rms - 10^(drc_threshold_dB/20))/drc_ratio;
    else
        compression_gain = signal_rms;
    end
    compressed_signal(i) = output_signal(i) / max(signal_rms, 1);
    
    % Update feedback buffer
    feedback_buffer = [compressed_signal(i), feedback_buffer(1:end-1)];
end

% Noise Reduction using Wiener Filter
if noise_reduction
    noise_estimate = 0.02; % Example noise level (adjust as needed)
    
    % Compute the power spectrum of the noisy signal
    N = length(compressed_signal);
    X = fft(compressed_signal);
    P = abs(X).^2 / N; % Power spectrum
    
    % Wiener filter: noise estimation and suppression
    noise_power = noise_estimate * ones(size(P)); % Constant noise estimate
    P_clean = max(P - noise_power, 0); % Subtract noise estimate and avoid negative values
    
    % Compute the inverse FFT to get the time-domain signal
    X_clean = sqrt(P_clean) .* exp(1i * angle(X));
    clean_signal = real(ifft(X_clean, 'symmetric'));
    
    % Ensure the clean signal has the same length as the input signal
    clean_signal = clean_signal(1:N);
else
    clean_signal = compressed_signal;
end

% Feedback Suppression (Dynamic)
% Instead of a fixed notch, we estimate feedback frequency dynamically
f0 = 1000; % Initial guess for frequency
Q = 30; % Quality factor

% Feedback suppression using an adaptive notch filter
% Adjust this to use a longer window for frequency estimation
windowSize = min(2*Fs, length(clean_signal)); % Use 2 seconds window or the length of the signal
for i = 1:length(clean_signal)-windowSize+1
    % Adaptive notch: use frequency estimation
    segment = clean_signal(i:i+windowSize-1);
    [~, freqIdx] = max(abs(fft(segment)));
    freq = (freqIdx - 1) * Fs / windowSize; % Convert index to frequency
    
    % Design notch filter
    notch_filter = designfilt('bandstopiir', 'FilterOrder', 2, ...
        'HalfPowerFrequency1', freq-50, 'HalfPowerFrequency2', freq+50, ...
        'DesignMethod', 'butter', 'SampleRate', Fs);
    
    % Apply the notch filter to the entire signal
    clean_signal = filtfilt(notch_filter, clean_signal);
end

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
plot(t, clean_signal);
title('Final Output Signal with Noise Reduction and Feedback Suppression');
xlabel('Time (sec)');
ylabel('Amplitude');

% Save the final output signal to a file
audiowrite('final_output_signal.wav', clean_signal, Fs);
disp('Simulation complete. The final output signal is saved as ''final_output_signal.wav''.');
