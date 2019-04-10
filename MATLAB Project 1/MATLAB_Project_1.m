%%% MATLAB Project 1

%%% Variables and Constants

F_1F = 1e9; % Intermediate Frequency
T_p = 0.5e-6; % End of time vector
B = 100e6; % Bandwidth
P = 1; % Power
R = 1000; % Distance
c = 3e8; % The Speed of Light
F_c = 8e9; % Center Frequency
alpha = 1; % Another amplitude
tau = (2 .* R) ./ c; % Time it takes to travel there and back
upconverting_freq = 7e9; % Frequency of the cosine

%%% Background

% Time array
t = 0:0.000000001: T_p;

% Chirp rate
gamma = B ./ T_p;

% Amplitude
 

% x_1F
x_1F = a .* cos((2 .* pi .* (F_1F - (0.5 .* B)) .* t) + (pi .* gamma .* t .* t));

% x(t) 
x = a .* cos(((2 .* pi) .* ((2 .* pi) .* (F_c - (0.5 .* B))) .* t) + (pi .* gamma .* t .* t));

% g(t)
g = alpha .* a .* cos((2 .* pi .* F_c .* t) - (2 .* pi * F_c .* tau) - (2 .* pi .* (-0.5) .* B .* t) + (2 .* pi .* 0.5 .* tau) + (pi .* gamma .* (t - tau) .* (t - tau)));

% Noise
% noise = sqrt(Pn ./ 2) .* (randn(L,1) + (1j .* randn(L,1)));

%%% Part 1

% Upconverting
upconvert = cos(2 .* pi .* upconfreq .* t);
x_1F_upconverted = upconvert .* x_1F;

% Bandpass filter

