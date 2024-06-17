close all 
%%% PART 1
% Define an array of stimulus intensities (z values)
z_values = linspace(-2, 1, 10); % Adjust the range and number of values as needed
IC = [0 0]; 
tspan = [0 200]; 
% Loop through each stimulus intensity
figure(1)
colorMap = jet(length(z_values));
for i = 1:length(z_values)
    % Set the current stimulus intensity
    z = z_values(i);

    currentColor = colorMap(i, :);

    % Solve the ODEs for the current stimulus intensity
    [t, xy] = ode15s(@(t, xy) fn(t, xy, z), tspan, IC);
    % Extract the solution for x
    x = xy(:, 1);
    % Plot the solution for x
    plot(t, x, 'DisplayName', sprintf('z = %.2f', z), 'Color', currentColor);
    hold on;
end

hold off;
grid on;
xlabel('Time');
ylabel('x');
legend('Location','BestOutside');
title('FitzHugh-Nagumo Model at Different Stimulus Intensities');

figure(2)
colorMap = jet(length(z_values));
for i = 1:length(z_values)
    % Set the current stimulus intensity
    z = z_values(i);

    currentColor = colorMap(i, :);

    % Solve the ODEs for the current stimulus intensity
    [t, xy] = ode15s(@(t, xy) fn(t, xy, z), tspan, IC);

    % Extract the solution for x
    y = xy(:, 2);
    % Plot the solution for x
    plot(t, y, 'DisplayName', sprintf('z = %.2f', z),'Color', currentColor);
    hold on;
end

hold off;
grid on;
xlabel('Time');
ylabel('y');
legend('Location','BestOutside');
title('FitzHugh-Nagumo Model at Different Stimulus Intensities');

figure(3)
% Use a color map for distinct colors
colorMap = winter(length(z_values));
for i = 1:length(z_values)
    % Set the current stimulus intensity
    z = z_values(i);

    % Get the color for the current iteration
    currentColor = colorMap(i, :);

    % Solve the ODEs for the current stimulus intensity
    [t, xy] = ode15s(@(t, xy) fn(t, xy, z), tspan, IC);
    
    % Extract the solution for x
    x = xy(:, 1);
    y = xy(:, 2);
    % Plot the solution for x
    plot(x, y, 'DisplayName', sprintf('z = %.2f', z), 'Color', currentColor);
    hold on;
    w = 1;
    a = 0.7;
    b = 0.8;
    x1 = linspace(-2.5, 2.5, 200); 
    y1 = (x1.^3)/3 - x1 - z; 
    plot(x1, y1, '--', 'DisplayName',sprintf('dx/dt=0, z = %.2f', z), 'Color', currentColor); 
    hold on; 
end

w = 1;
a = 0.7;
b = 0.8;
z = 0; 
x2 = linspace(-2.5, 2.5, 200); 
y2 = a/b - (w^2)*x2./b ; 
plot(x2,y2, 'DisplayName','dy/dt = 0'); 
hold on; 

hold off;
grid on;
xlabel('x');
ylabel('y');
legend('Location','BestOutside');
title('FitzHugh-Nagumo Phase-Space Plot at Different Stimulus Intensities');

figure(4)
t_span = [0, 200];
y0 = [0, 0];
z_values = linspace(-1.5, 0, 1000); % can make it larger for more fine results
intrinsic_freq = []; % empty array to store the intrinsic frequencies extracted from each z value
frequencies = 0.01:0.001:1; % create frequencies to look at. Will expect something in between 0 and 1 Hz for neurons (tried 2, got similar results so condensed it)
scales = 1 ./ frequencies;
fs = 1000; % convert the frequency to scale for pywavelets function used later

%wavelet = cmorwavf(-1, 1, 1000, 1, 1);

for i = z_values
    % Solve the ODEs for the current stimulus intensity
    [t, xy] = ode15s(@(t, xy) fn(t, xy, i), tspan, IC);
    % Extract the solution for x
    x = xy(:, 1);
    % utilize complex morlet wavelet transform for x solution and set scale
    [coeffs, freqs] = cwt(x, scales, 'cmor', fs); % retrieve the coefficients and the respective frequencies from the cwt
    [~, idx] = max(sum(abs(coeffs), 2)); % find the index where the sum of the absolute value of the coefficient in that row (frequency scale) is largest
    dom_freq = freqs(idx); % find the frequency values corresponding to the previous sum; this is the dominant or intrinsic freq
    intrinsic_freq = [intrinsic_freq, dom_freq]; % add that freq to the empty list
end

% Visualize and figure out exact zeros for determining k1 and k2 range
figure;
plot(z_values, intrinsic_freq);
xlabel('Stimulus Intensity (Z)');
ylabel('Frequency');
title('Intrinsic Frequency v. Stimulus Intensity');

%%% PART 2


% FitzHugh-Nagumo ODE function with additional parameter z
function dxydt = fn(t, xy, z)
    alpha = 3;
    w = 1;
    a = 0.7;
    b = 0.8;
    x = xy(1);
    y = xy(2);
    dxydt = [alpha * (y + x - (x^3)/3 + z);
             (-1/alpha) * ((w^2) * x - a + b * y)];
end


function coupled_oscillator = f_coupled(t, xy, k1, k2, c)
    alpha = 3;
    w = 1;
    a = 0.7;
    b = 0.8;
    x1 = xy(1);
    y1 = xy(2);
    x2 = xy(3);
    y2 = xy(4);
    coupled_oscillator = [alpha * (y1 + x1 - (x1.^3)/3 + k1 + c*x2);
             (-1/alpha) * ((w^2) * x1 - a + b * y1);
             alpha * (y2 + x2 - (x2.^3)/3 + k2 + c*x1);
             (-1/alpha) * ((w^2) * x2 - a + b * y2)]; 
end 