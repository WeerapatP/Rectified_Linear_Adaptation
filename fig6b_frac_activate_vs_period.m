set(0,'defaultAxesFontName', 'Arial')
set(0,'defaultTextFontName', 'Arial')
num_traces = 12;

activation_threshold = 0.3;
max_signal = 0.1;
offset = max_signal*5/10;
amplitude = max_signal*4/10;
noise_str = 0;
flag_plot = false;
dt = 0.0002;
final_time = 150;
period = logspace(-1, 2, num_traces);
noise_corr_time = period;
ensemblesize = 50;
sigma = 0.01;
times = 0:dt:final_time;
time_activate_avg = zeros(num_traces, 1);
time_peak     = zeros(ensemblesize, num_traces);
frac_activate = zeros(num_traces);
flag_rectified = true;
flag_square = false;

for i = 1:num_traces
    [~, ~, ~, ~, time_activate_avg(i, 1), time_peak(:,i),  frac_activate(i)] = simulate_abstract_NFkB(offset, amplitude, period(i), noise_str, noise_corr_time(i), flag_plot, flag_rectified, flag_square, dt, final_time, ensemblesize, sigma, activation_threshold); 
end

%%
figure()
semilogx(period', time_activate_avg, '-o', 'LineWidth', 2)
xlabel('Period (\tau_0)')
ylabel('Average Activation Time (min)')
set(findall(gcf,'-property','FontSize'),'FontSize',20)
set(gca,'box','off')

figure()
semilogx(period, frac_activate, '-o', 'LineWidth', 2)
xlabel('Period (\tau_0)')
ylabel("Fraction Activated (NFkB > " + num2str(activation_threshold) + ")")
set(findall(gcf,'-property','FontSize'),'FontSize',20)
set(gca,'box','off')

figure()
errorbar(log10(period), mean(time_peak), sqrt(var(time_peak)), 'LineWidth', 2)
xlabel('Period (\tau_0)')
ylabel('Peak Time (min)')
xticks([-1, -0.5, 0])
xticklabels({'10^{-1}', '10^{-0.5}', '10^0'})
set(findall(gcf,'-property','FontSize'),'FontSize',20)
set(gca,'box','off')