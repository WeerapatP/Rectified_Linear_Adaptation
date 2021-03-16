
set(0,'defaultAxesFontName', 'Arial')
set(0,'defaultTextFontName', 'Arial')
num_traces = 2;

activation_threshold = 0.25;
max_signal = 0.1;
noise_amp = max_signal/4;
offset = [max_signal, noise_amp];
amplitude = 0;
noise_str = [0, noise_amp];
flag_plot = false;
dt = 0.0002;
final_time = 200;
period = 4/60*ones(1, 2);
noise_corr_time = period;
ensemblesize = 10;
sigma = 0.01;
flag_rectified = true;
flag_square = false;

times = 0:dt:final_time;
trace_inputs = zeros(num_traces, length(times));
trace_ys     = zeros(num_traces, length(times));
trace_IKKas  = zeros(num_traces, length(times));
trace_NFkBs  = zeros(num_traces, length(times), ensemblesize);
time_activate_avgs = zeros(num_traces, 1);
for i = 1:num_traces
    [trace_input, trace_y, trace_IKKa, trace_NFkB, time_activate_avg, ~] = simulate_abstract_NFkB(offset(i), amplitude, period(i), noise_str(i), noise_corr_time(i), flag_plot, flag_rectified, flag_square, dt, final_time, ensemblesize, sigma, activation_threshold);
    trace_inputs(i, :) = trace_input;
    trace_ys(i, :) = trace_y;
    trace_IKKas(i, :) = trace_IKKa;
    trace_NFkBs(i, :, :) = trace_NFkB; 
    time_activate_avgs(i) = time_activate_avg;
end

%%
time_mul_factor = 120;
set_y_max = 0.17;
set_IKKa_max = 0.35;
set_NFkB_max = 0.65;
ylabel_texts = ["Noisy", "Noiseless"];
%%
figure('rend','painters','pos',[0 0 500 1200])
subplot(3, 1, 1)
hold all;
plot(times, trace_inputs(2, :), 'Color', 'Red', 'LineWidth', 2)
plot(times, trace_inputs(1, :), 'Color', 'Blue', 'LineWidth', 2)
hold off;
xlim([0, noise_corr_time(i)*time_mul_factor])
% ylim([0, max_signal])
xlabel('time (\tau_0)')
ylabel('TNF-\alpha (Input)')

subplot(3, 1, 2)
hold all
for j = 1:ensemblesize
    plot(times, reshape(trace_NFkBs(2, :, j), 1, []),'Color', 'Red', 'LineWidth', 1)
end
ylim([0, set_NFkB_max])
xlim([0, final_time])
hold off
title('Noisy Input')
ylabel('NF-kB')
xlabel('time (min)')
    
subplot(3, 1, 3)
hold all
for j = 1:ensemblesize
    plot(times, reshape(trace_NFkBs(1, :, j), 1, []),'Color', 'Blue', 'LineWidth', 1)
end
ylim([0, set_NFkB_max])
xlim([0, final_time])
hold off
title('Noiseless Input')
ylabel('NF-kB')
xlabel('time (min)')
set(findall(gcf,'-property','FontSize'),'FontSize', 22)
