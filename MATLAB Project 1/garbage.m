T_p =0.5e-6;
min = -0.5 + (0.5 * T_p);
max = 0.5 + (0.5 * T_p);
a = T_p .* rectangularPulse(min,max,t);

plot(t,a)
ylim([0 (2 * T_p)])