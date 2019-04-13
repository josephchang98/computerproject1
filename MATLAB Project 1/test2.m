clear all;
close all;

%%% MATLAB Project 1

%%% Variables and Constants

F_1F = 1e9; % Intermediate Frequency
T_p = 0.5e-6; % End of time vector
B = 100e6; % Bandwidth
P = 1; % Power
A = sqrt(2); % Amplitude = sqrt(2*Power)
R = 1000; % Distance
c = 3e8; % The Speed of Light
F_c = 8e9; % Center Frequency
alpha = 1; % Another amplitude
tau = (2 .* R) ./ c; % Time it takes to travel there and back
upconverting_freq = 7e9; % Frequency of the cosine
sflag = 20e9; % Sampling
Pn = 10; % Noise power
sampling_rate = 20e9;
F_s = 1 ./ sampling_rate; % Period
n = 2; % Downsampling -- CHANGE THIS IT IS WRONG

%%% Background

% Time array
t = 0:F_s:T_p;
% Chirp rate
gamma = B ./ T_p;
% Amplitude
a = rectangularPulse(0,T_p,t);
% x_1F
x_1F = A .* a .* cos((2 .* pi .* (F_1F - (0.5 .* B)) .* t) + (pi .* gamma .* t .* t));
% x(t) 
x = a .* cos(((2 .* pi) .* ((2 .* pi) .* (F_c - (0.5 .* B))) .* t) + (pi .* gamma .* t .* t));
% g(t)
g = alpha .* a .* cos((2 .* pi .* F_c .* t) - (2 .* pi * F_c .* tau) - (2 .* pi .* (-0.5) .* B .* t) + (2 .* pi .* 0.5 .* tau) + (pi .* gamma .* (t - tau) .* (t - tau)));

%%% Part 1
% Upconverting
upconvert = cos(2 .* pi .* upconverting_freq .* t);
x_1F_upconverted = upconvert .* x_1F;
F = -1e14:20e9:1e14;
F_f = linspace(-1e10, 1e10, length(F));
figure(1)
plot(F_f,fftshift(abs(fft(x_1F))))
xlabel('Frequency (Hz)');
ylabel('Signal Strength (V/m)');

figure(2)
plot(F_f,fftshift(abs(fft(x_1F_upconverted))))
xlabel('Frequency (Hz)');
ylabel('Signal Strength (V/m)');

% Bandpass filter
T_1 = -250e-9; % Not sure about this
T_2 = 250e-9; % Not sure about this
t_BPF = T_1:F_s:T_2;
W = 4 .* B;
w = W .* transpose(hamming(length(t_BPF)));
h = w .* sinc(W .* t_BPF);
h_BPF = h .* cos(2 .* pi .* F_c .* t_BPF);
figure(3)
plot(F_f,fftshift(abs(fft(h_BPF))))
xlabel('Frequency (Hz)');
ylabel('Signal Strength (V/m)');

% Convolution
x_BPF = conv(h_BPF,x_1F_upconverted);
figure(4)
F_conv = -2e14:20e9:2e14;
F_f_conv = linspace(-2e10, 2e10, length(F_conv));
figure(4)
plot(F_f_conv, fftshift(abs(x_BPF)))
xlabel('Frequency (Hz)');
ylabel('Signal Strength (V/m)');

t_conv = 0:F_s:2*T_p; % time array for x_BPF; if you plot it, it is centered at 0.5...delay or is t_conv wrong?

% Scaling and Truncation -- a problem for future allyson lol
x_BPF = (sqrt(2 .* P)./max(x_BPF)) .* x_BPF; % Very good! Thank you Joe

% New Time Arrays
delta_R = c / (2 .* B);
R_min = R - (0.5 .* delta_R);
R_max = R + (0.5 .* delta_R);
new_min = (2 .* R_min) ./ c;
new_max = ((2 .* R_max) ./ c) + T_p + T_2;
t_new = new_min:F_s:new_max;

% Signal Interpolation
x_interp = interp1(t_conv + tau, x_BPF, t_new, 'linear', 'extrap'); % i assume I'm shifting t_conv by tau here

% Downconversion
cos_new = cos(2 .* pi .* F_c .* t_new);
sin_new = sin(2 .* pi .* F_c .* t_new);
x_downconverted_cos = x_interp .* cos_new;
x_downconverted_sin = x_interp .* sin_new;
figure(5)
subplot(2,1,1);
plot(t_new,x_downconverted_cos);
subplot(2,1,2);
xlabel('Time (Seconds)');
ylabel('Signal Strength (V/m)');
plot(t_new,x_downconverted_sin);
xlabel('Time (Seconds)');
ylabel('Signal Strength (V/m)');

% Low Pass Filter
W_l = 2 .* (8e9); % want to keep 8GHz freq
w_l = W .* transpose(hamming(length(t_BPF)));
h_LBP =  w .* sinc(W .* t_BPF); % not sure about t_BPF
h = h_LBP;
D1_conv = conv(h_LBP, x_downconverted_cos);
D2_conv = conv(h_LBP, x_downconverted_sin);
figure(6)
subplot(2,1,1);
%plot(F_f_conv,fftshift(abs(fft(D1_conv))));
xlabel('Frequency (Hz)');
ylabel('Signal Strength (V/m)');
subplot(2,1,2);
%plot(F_f_conv,fftshift(abs(fft(D2_conv))));
xlabel('Frequency (Hz)');
ylabel('Signal Strength (V/m)');

% Scaling and Truncation-- YET ANOTHER problem for future allyson lol
x_downconverted_cos_lowpass = (sqrt(2 .* P) ./ max(x_downconverted_cos)) .* x_downconverted_cos;
x_downconverted_sin_lowpass = (sqrt(2 .* P) ./ max(x_downconverted_sin)) .* x_downconverted_sin;

% Combination
x_downconverted_lowpass = x_downconverted_cos_lowpass + x_downconverted_sin_lowpass;
figure(7)
%plot(F_f_conv, fftshift(abs(fft(x_downconverted_lowpass))));
xlabel('Frequency (Hz)');
ylabel('Signal Strength (V/m)');

% Downsampling
x_downconverted_lowpass_C = downsample(x_downconverted_lowpass,n);

% Adding Noise
L = length(t_new);
noise = sqrt(Pn ./ 2) .* (randn(L,1) + (1j .* randn(L,1)));
x_noise = x_downconverted_lowpass_C + noise; 
figure(8)
%plot(t_new,x_noise);
%plot(F_f,fftshift(abs(fft(x_noise))));
%plot(F_f_conv,x_noise)
xlabel('Frequency (Hz)');
ylabel('Signal Strength (V/m)');

% Baseband copy of waveform for matched filter
x_ = a .* exp(1j .* 2 .* pi .* (0.5 .* B) .* t + 1j .* pi .* gamma .* t .* t);
x_S = conj(fliplr(x_));
figure(9)
plot(F_f,fftshift(abs(fft(x_S))));
xlabel('Frequency (Hz)');
ylabel('Signal Strength (V/m)');
