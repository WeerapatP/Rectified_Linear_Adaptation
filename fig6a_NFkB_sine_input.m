set(0,'defaultAxesFontName', 'Arial')
set(0,'defaultTextFontName', 'Arial')
num_traces = 3;
num_col = 3;

activation_threshold = 0.3;
max_signal = 0.1;
offset = max_signal*5/10;
amplitude = max_signal*4/10;
noise_str = 0;
flag_plot = false;
dt = 0.0002;
final_time = 150;
period = 2.^linspace(4, 9, num_traces)*1/60;
noise_corr_time = period;
ensemblesize = 10;
sigma = 0.01;
flag_square = false;
flag_rectified = true;

times = 0:dt:final_time;
trace_inputs = zeros(num_traces, length(times));
trace_ys     = zeros(num_traces, length(times));
trace_IKKas  = zeros(num_traces, length(times));
trace_NFkBs  = zeros(num_traces, length(times), ensemblesize);
time_activate_avgs = zeros(num_traces, 1);
frac_activates     = zeros(num_traces, 1);

for i = 1:num_traces
    [trace_input, trace_y, trace_IKKa, trace_NFkB, time_activate_avg, frac_activate] = simulate_abstract_NFkB(offset, amplitude, period(i), noise_str, noise_corr_time(i), flag_plot, flag_rectified, flag_square, dt, final_time, ensemblesize, sigma, activation_threshold);
    trace_inputs(i, :) = trace_input;
    trace_ys(i, :) = trace_y;
    trace_IKKas(i, :) = trace_IKKa;
    trace_NFkBs(i, :, :) = trace_NFkB; 
    time_activate_avgs(i) = time_activate_avg;
%     frac_activates(i) = frac_activate;
end

%%
max_time = final_time;
time_mul_factor = 1;
max_y    = 1.1*max(max(trace_ys));
max_IKKa = 1.1*max(max(trace_IKKas));
max_NFkB = 1.1*max(max(max(trace_NFkBs)));

figure('rend','painters','pos',[0 0 1500 1000])
for i = 1:num_traces
    subplot(num_traces, num_col, i)
    plot(times, trace_inputs(i, :), 'Color', [87, 108, 67]/255, 'LineWidth', 2)
%     title({"Period "; num2str(round(period(i), 2)) + "\tau_0"})
    title(num2str(round(period(i), 2)) + "\tau_0")
    xlim([0, max(period)*time_mul_factor])
    xlabel('Time (\tau_0)')
    if i == 1
        ylabel('Input')
    end
    set(gca,'box','off')
    
    subplot(num_traces, num_col, i + 3)
    plot(times, trace_IKKas(i, :), 'Color', [87, 108, 67]/255, 'LineWidth', 2)
    ylim([0, max_IKKa])
    xlim([0, max_time])
    xlabel('Time (min)')
    if i == 1
        ylabel('IKKa')
    end
    set(gca,'box','off')
    
    subplot(num_traces, num_col, i + 6)
    hold all
    for j = 1:ensemblesize
        plot(times, reshape(trace_NFkBs(i, :, j), 1, []),'Color', [87, 108, 67]/255, 'LineWidth', 1)
        ylim([0, max_NFkB])
    end
    ylim([0, max_NFkB])
    xlim([0, max_time])
    hold off
    xlabel('Time (min)')
    if i == 1
        ylabel('NF-kB')
    end
    set(gca,'box','off')
end
set(findall(gcf,'-property','FontSize'),'FontSize',20)