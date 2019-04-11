clear all;
close all;

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
sflag = 20e9; % Sampling
Pn = 10; % Noise power
sampling_rate = 20e9;
F_s = 1 ./ sampling_rate; % Sampling rate

%%% Background

% Time array
t = 0:F_s: T_p;

% Chirp rate
gamma = B ./ T_p;

% Amplitude
a = rectangularPulse(0,T_p,t);

% x_1F
x_1F = a .* cos((2 .* pi .* (F_1F - (0.5 .* B)) .* t) + (pi .* gamma .* t .* t));

% x(t) 
x = a .* cos(((2 .* pi) .* ((2 .* pi) .* (F_c - (0.5 .* B))) .* t) + (pi .* gamma .* t .* t));

% g(t)
g = alpha .* a .* cos((2 .* pi .* F_c .* t) - (2 .* pi * F_c .* tau) - (2 .* pi .* (-0.5) .* B .* t) + (2 .* pi .* 0.5 .* tau) + (pi .* gamma .* (t - tau) .* (t - tau)));

% Noise
% L = idk
% noise = sqrt(Pn ./ 2) .* (randn(L,1) + (1j .* randn(L,1)));

%%% Part 1

% Upconverting
upconvert = cos(2 .* pi .* upconverting_freq .* t);
x_1F_upconverted = upconvert .* x_1F;
F = -1e14:20e9:1e14;
figure(1)
plot(F,fftshift(abs(fft(x_1F_upconverted))))

% Bandpass filter
T_1 = -250e-9; % Not sure about this
T_2 = 250e-9; % Not sure about this
t_BPF = T_1:F_s:T_2;
W = 4 .* B;
w = W .* transpose(hamming(length(t_BPF)));
h = w .* sinc(W .* t_BPF);
h_BPF = h .* cos(2 .* pi .* F_c .* t_BPF);
figure(2)
plot(F,fftshift(abs(fft(h_BPF))))

% Convolution
x_BPF = conv(h_BPF,x_1F_upconverted);

% Scaling and Truncation -- a problem for future allyson lmao


% New Time Arrays
delta_R = c / (2 .* B);
R_min = R - (0.5 .* delta_R);
R_max = R + (0.5 .* delta_R);
new_min = (2 .* R_min) ./ c;
new_max = ((2 .* R_max) ./ c) + T_p + T_2;
t_new = new_min:F_s:new_max;

% Signal Interpolation

% Downconversion
cos_new = cos(2 .* pi .* F_c .* t_new);
sin_new = sin(2 .* pi .* F_c .* t_new);

% Scaling and Truncation-- YET ANOTHER problem for future allyson lmao
