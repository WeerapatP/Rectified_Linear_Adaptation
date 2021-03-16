activation_threshold = 0.3;
max_signal = 0.1;
offset = 0;
noise_str = 0;
amplitude = max_signal;
flag_plot = true;
dt = 0.0002;
final_time = 150;
period = 60/60;
noise_corr_time = period;
ensemblesize = 5;
sigma = 0.0;
flag_square = true;

flag_rectified = true;
simulate_abstract_NFkB(offset, amplitude, period, noise_str, noise_corr_time, flag_plot, flag_rectified, flag_square, dt, final_time, ensemblesize, sigma, activation_threshold);

flag_rectified = false;
simulate_abstract_NFkB(offset, amplitude, period, noise_str, noise_corr_time, flag_plot, flag_rectified, flag_square, dt, final_time, ensemblesize, sigma, activation_threshold);