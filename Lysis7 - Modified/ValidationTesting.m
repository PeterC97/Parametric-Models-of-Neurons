% Define parameters
Fs = 1000;          % Sampling frequency (Hz)
t = 0:1/Fs:1;       % Time vector (1 second duration)
A = 1;              % Amplitude

% Generate sinusoidal signals
signal_6Hz = A * sin(2 * pi * 6 * t);
signal_10Hz = A * sin(2 * pi * 10 * t);
signal_20Hz = A * sin(2 * pi * 20 * t);

% Superimpose the signals
combined_signal = signal_6Hz + signal_10Hz + signal_20Hz;

% Plot the combined signal
figure;
subplot(2,1,1);
plot(t, combined_signal);
title('Combined Sinusoidal Signals');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% Perform FFT
N = length(combined_signal);                % Length of the signal
frequencies = (0:N-1)*(Fs/N);               % Frequency vector
fft_signal = fft(combined_signal)/N;        % Compute FFT
fft_signal_abs = abs(fft_signal(1:N/2+1));  % Take the absolute value of the FFT result

% Plot the FFT result
subplot(2,1,2);
plot(frequencies(1:N/2+1), fft_signal_abs);
title('FFT of Combined Signal');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;


alphas = linspace(0, 1, 21);  % Include 1 and exclude 0
alphas = alphas(2:end); 
Ls = linspace(1,10,10);
num_best = 5; % Number of best indexes you want to find
best_values = Inf(1, num_best); % Initialize with Inf
best_alphas = zeros(1, num_best);
best_Ls = zeros(1, num_best);

for i = 1:length(alphas)
    alpha = alphas(i);
    for j = 1:length(Ls)
        L = Ls(j);
        [Cest, Kest, Pred, NMSE] = LET_1(low_noise, x1l, alpha, L, Q, Nfig);
        sum_NMSE = sum(abs(NMSE));
        [max_value, max_index] = max(best_values);
        if sum_NMSE < max_value
            best_values(max_index) = sum_NMSE;
            best_alphas(max_index) = i;
            best_Ls(max_index) = j;
        end 
    end
end

% Output the lowest 5 alphas and Ls
for k = 1:num_best
    disp(['Alpha ', num2str(k), ': ', num2str(alphas(best_alphas(k))), ', L ', num2str(k), ': ', num2str(Ls(best_Ls(k)))]);
end