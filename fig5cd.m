
set(0,'defaultAxesFontName', 'Arial')
set(0,'defaultTextFontName', 'Arial')
num_traces = 2;

activation_threshold = 0.3;
max_signal = 0.1;
offset = max_signal*5/10;
amplitude = 0;
noise_str = [max_signal*4/10, 0];
flag_plot = false;
dt = 0.0002;
final_time = 150;
period = 4/60*ones(1, 2);
noise_corr_time = period;
ensemblesize = 10;
sigma = 0.01;

times = 0:dt:final_time;
trace_inputs = zeros(num_traces, length(times));
trace_ys     = zeros(num_traces, length(times));
trace_IKKas  = zeros(num_traces, length(times));
trace_NFkBs  = zeros(num_traces, length(times), ensemblesize);
time_activate_avgs = zeros(num_traces, 1);
frac_activates     = zeros(num_traces, 1);
flag_square = false;

% fig3c
flag_rectified = false;

for i = 1:num_traces
    [trace_input, trace_y, trace_IKKa, trace_NFkB, time_activate_avg, frac_activate] = simulate_abstract_NFkB(offset, amplitude, period(i), noise_str(i), noise_corr_time(i), flag_plot, flag_rectified, flag_square, dt, final_time, ensemblesize, sigma, activation_threshold);
    trace_inputs(i, :) = trace_input;
    trace_ys(i, :) = trace_y;
    trace_IKKas(i, :) = trace_IKKa;
    trace_NFkBs(i, :, :) = trace_NFkB; 
    time_activate_avgs(i) = time_activate_avg;
%     frac_activates(i) = frac_activate;
end

%%
time_mul_factor = 120;
set_y_max = 0.17;
set_IKKa_max = 0.35;
set_NFkB_max = 0.65;
ylabel_texts = ["Noisy Signal", "Constant Signal"];

figure('rend','painters','pos',[0 0 1500 600])
colors=["Red", "Blue"];
for i = 1:num_traces
    subplot(num_traces, 4, (i-1)*4 + 1)
    plot(times, trace_inputs(i, :), 'Color', char(colors(i)), 'LineWidth', 2)
    ylabel(ylabel_texts(i))
    xlim([0, noise_corr_time(i)*time_mul_factor])
    ylim([0, max_signal])
    if i == 1
        title('Input')
    elseif i == num_traces
        xlabel('time (\tau_0)')
    end
    set(gca,'box','off')
    
    subplot(num_traces, 4, (i-1)*4 + 2)
    plot(times, trace_ys(i, :), 'Color', char(colors(i)), 'LineWidth', 2)
    ylim([0, set_y_max])
    xlim([0, noise_corr_time(i)*time_mul_factor])
    if i == 1
        title('Y')
    elseif i == num_traces
        xlabel('time (\tau_0)')
    end
    set(gca,'box','off')
    
    subplot(num_traces, 4, (i-1)*4 + 3)
    plot(times, trace_IKKas(i, :), 'Color', char(colors(i)), 'LineWidth', 2)
    ylim([0, set_IKKa_max])
    xlim([0, final_time])
    if i == 1
        title('IKKa')
    elseif i == num_traces
        xlabel('time (min)')
    end
    set(gca,'box','off')
    
    subplot(num_traces, 4, (i-1)*4 + 4)
    hold all
    for j = 1:ensemblesize
        plot(times, reshape(trace_NFkBs(i, :, j), 1, []),'Color', char(colors(i)), 'LineWidth', 1)
    end
    ylim([0, set_NFkB_max])
    xlim([0, final_time])
    hold off
    if i == 1
        title('NF-kB')
    elseif i == num_traces
        xlabel('time (min)')
    end
    set(gca,'box','off')
end
set(findall(gcf,'-property','FontSize'),'FontSize', 22)

%% fig3d
flag_rectified = true;

for i = 1:num_traces
    [trace_input, trace_y, trace_IKKa, trace_NFkB, time_activate_avg, frac_activate] = simulate_abstract_NFkB(offset, amplitude, period(i), noise_str(i), noise_corr_time(i), flag_plot, flag_rectified, flag_square, dt, final_time, ensemblesize, sigma, activation_threshold);
    trace_inputs(i, :) = trace_input;
    trace_ys(i, :) = trace_y;
    trace_IKKas(i, :) = trace_IKKa;
    trace_NFkBs(i, :, :) = trace_NFkB; 
    time_activate_avgs(i) = time_activate_avg;
end

%%
time_mul_factor = 120;
max_y    = 1.1*max(max(trace_ys));
max_IKKa = 1.1*max(max(trace_IKKas));
max_NFkB = 1.1*max(max(max(trace_NFkBs)));
set_y_max = 0.17;
set_IKKa_max = 0.35;
set_NFkB_max = 0.65;
ylabel_texts = ["Noisy Signal", "Constant Signal"];

figure('rend','painters','pos',[0 0 1500 600])
colors=["Red", "Blue"];
for i = 1:num_traces
    subplot(num_traces, 4, (i-1)*4 + 1)
    plot(times, trace_inputs(i, :), 'Color', char(colors(i)), 'LineWidth', 2)
    ylabel(ylabel_texts(i))
    xlim([0, noise_corr_time(i)*time_mul_factor])
    ylim([0, max_signal])
    if i == 1
        title('Input')
    elseif i == num_traces
        xlabel('time (\tau_0)')
    end
    set(gca,'box','off')
    
    subplot(num_traces, 4, (i-1)*4 + 2)
    plot(times, trace_ys(i, :), 'Color', char(colors(i)), 'LineWidth', 2)
    ylim([0, set_y_max])
    xlim([0, noise_corr_time(i)*time_mul_factor])
    if i == 1
        title('Y')
    elseif i == num_traces
        xlabel('time (\tau_0)')
    end
    set(gca,'box','off')
    
    subplot(num_traces, 4, (i-1)*4 + 3)
    plot(times, trace_IKKas(i, :), 'Color', char(colors(i)), 'LineWidth', 2)
    ylim([0, set_IKKa_max])
    xlim([0, final_time])
    if i == 1
        title('IKKa')
    elseif i == num_traces
        xlabel('time (min)')
    end
    set(gca,'box','off')
    
    subplot(num_traces, 4, (i-1)*4 + 4)
    hold all
    for j = 1:ensemblesize
        plot(times, reshape(trace_NFkBs(i, :, j), 1, []),'Color', char(colors(i)), 'LineWidth', 1)
    end
    ylim([0, set_NFkB_max])
    xlim([0, final_time])
    hold off
    if i == 1
        title('NF-kB')
    elseif i == num_traces
        xlabel('time (min)')
    end
    set(gca,'box','off')
end
set(findall(gcf,'-property','FontSize'),'FontSize', 22)