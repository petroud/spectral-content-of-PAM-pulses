clear all;
close all;
clc;

%----------------------------%
% Author: Petrou Dimitrios
% Year: 2023  
% TU of Crete
% Telecommunication Systems I
%----------------------------%


                         %---Part B---%
%----------------------------------------------------------------%
%---B1---%
% Parameters
F0 = 1;  % Frequency
t = linspace(0, 1, 1000);  % Time vector

% Generate 5 random samples
numSamples = 5;
X = randn(numSamples, 1);  % Normally distributed random numbers
Phi = 2*pi*rand(numSamples, 1);  % Uniformly distributed random numbers between 0 and 2ð

% Evaluate the equation for each sample
Y = X .* cos(2*pi*F0*t + Phi);

% Plot the results
figure;
hold on;

for i = 1:numSamples
    plot(t, Y(i, :));
end

xlabel('Time');
ylabel('Y(t)');
title('Multiple Implementations of Y(t) = X cos(2ðF0t + Ö)');
legend('Sample 1', 'Sample 2', 'Sample 3', 'Sample 4', 'Sample 5');

hold off;

