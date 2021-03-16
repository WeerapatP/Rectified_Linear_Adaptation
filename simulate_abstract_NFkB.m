function [inputs, trace_y, trace_IKKa, trace_NFkB, time_activate_avg, time_peak, frac_activate] = simulate_abstract_NFkB(offset, amplitude, period, noise_str, noise_corr_time, flag_plot, flag_rectified, flag_square, dt, final_time, ensemblesize, sigma, activation_threshold)
    % Author: Weerapat Pittayakanchit
    set(0,'defaultAxesFontName', 'Arial')
    set(0,'defaultTextFontName', 'Arial')
    %% Set the simulation specfic details
    t_initial = 0;
    times = 0:dt:final_time;
    inputs = offset*ones(length(times), 1);
    noise = 0;
    t_next = 0;
    
    for index = 1:length(times)
        time  = times(index);
        if (noise_str ~= 0) && (time >= t_next)
            t_next = t_next + exprnd(noise_corr_time);
            noise  = noise_str*(2*rand() - 1);
        end
            
        if time > t_initial
            % square wave
            if flag_square
                inputs(index) = max(offset + amplitude*(1 + square(time*(2*pi)/period))/2 + noise, 0);
            else
                % cosine input
                inputs(index) = max(offset + amplitude*sin(2*pi*time/period) + noise, 0);
            end
        else
            inputs(index) = max(offset + noise, 0);
        end
    end
    
    %% Set the initial conditions of the variables in the model
    Nn   = 0.0903*ones(ensemblesize, 1);
    Im   = 0.4941*ones(ensemblesize, 1);
    I    = 3.5348*ones(ensemblesize, 1);
    IKKa = zeros(ensemblesize, 1);
    y    = zeros(ensemblesize, 1);
    x    = zeros(ensemblesize, 1);
    time_activate = zeros(ensemblesize, 1);
    
    %% Set observable variables to be plotted
    trace_NFkB = zeros(length(times), ensemblesize);
    trace_IKKa = zeros(length(times), 1);
    trace_Im   = zeros(length(times), 1);
    trace_I    = zeros(length(times), 1);
    trace_y    = zeros(length(times), 1);
    trace_x    = zeros(length(times), 1);
    
    %% Parameters of the models
    kNin = 5.4;         klin = 0.018;   kt = 1.03;
    ktl  = 0.24;        KI  = 0.035;    KN = 0.029; 
    gamma_m  = 0.017;   alpha = 1.05;   Ntot = 1;
    alpha_x = 10;
    alpha_IKKa = 0.1;
    k_y = 30;
    beta_x  = 60;
    beta_TNF = 60.6;
    beta_y = 60;
    y0 = 0.001;
    if flag_rectified
        c0 = 0;
    else
        c0 = 0.1;
    end
    n=10;
        
    %% Initial condition
    pre_condition_times = 0:dt:final_time;
    for index = 1:length(pre_condition_times)
        time = pre_condition_times(index);
        TNF  = 0;
        
        dNn     = dt*(  kNin * (Ntot - Nn) * KI./(KI + I)    - klin * I .* Nn ./(Nn + KN));
        dIm     = dt*(  kt * Nn.^2                           - gamma_m * Im);
        dI      = dt*(  ktl * Im                             - alpha * IKKa .* (Ntot - Nn) .* I ./(KI + I));
        
        hill = ((y.^n)./(y0.^n + y.^n));
        dIKKa   = -dt*alpha_IKKa*(IKKa - k_y*(y-c0));
        dy      =  dt*(-beta_y*(y-c0) - beta_x*x.*hill + beta_TNF*TNF);
        dx      = -dt*alpha_x*(x - TNF);
        
        Nn      = Nn + dNn;
        Im      = Im + dIm;
        I       = I  + dI;
        IKKa    = IKKa + dIKKa;
        if sigma ~= 0
            Nn = Nn + normrnd(0.0,sigma, ensemblesize, 1)*sqrt(dt);
            Im = Im + normrnd(0.0,sigma, ensemblesize, 1)*sqrt(dt);
            I  =  I + normrnd(0.0,sigma, ensemblesize, 1)*sqrt(dt);
        end
        if flag_rectified
            y       = poslin(y + dy);
        else
            y       = y + dy;
        end
        x       = x + dx;
    end
    
    
    %% Simulate
    for index = 1:length(times)
        time  = times(index);
        input = inputs(index);
        TNF   = input;
        
        dNn     = dt*(  kNin * (Ntot - Nn) * KI./(KI + I)    - klin * I .* Nn ./(Nn + KN));
        dIm     = dt*(  kt * Nn.^2                           - gamma_m * Im);
        dI      = dt*(  ktl * Im                             - alpha * IKKa .* (Ntot - Nn) .* I ./(KI + I));
        
        hill = ((y.^n)./(y0.^n + y.^n));
        dIKKa   = -dt*alpha_IKKa*(IKKa - k_y*(y-c0));
        dy      =  dt*(-beta_y*(y-c0) - beta_x*x.*hill + beta_TNF*TNF);
        dx      = -dt*alpha_x*(x - TNF);
        
        Nn      = Nn + dNn;
        Im      = Im + dIm;
        I       = I  + dI;
        IKKa    = IKKa + dIKKa;
        if sigma ~= 0
            Nn = Nn + normrnd(0.0,sigma, ensemblesize, 1)*sqrt(dt);
            Im = Im + normrnd(0.0,sigma, ensemblesize, 1)*sqrt(dt);
            I = I + normrnd(0.0,sigma, ensemblesize, 1)*sqrt(dt);
        end
        if flag_rectified
            y       = poslin(y + dy);
        else
            y       = y + dy;
        end
        x       = x + dx;
      
        trace_NFkB(index, :) = Nn;
        trace_Im(index)   = Im(1);
        trace_I(index)    = I(1);
        trace_IKKa(index) = IKKa(1);
        trace_y(index)    = y(1);
        trace_x(index)    = x(1);
    end
    
    %% Plot
    set(0,'defaultAxesFontName', 'Arial')
    set(0,'defaultTextFontName', 'Arial')
%     if flag_plot == true
%         figure()
% 
%         subplot(3,2,1)
%         plot(times, inputs)
%         ylabel('Input')
%         ylim([0, 1.5*max(inputs)])
%         
% 
%         subplot(3,2,2)
%         plot(times, trace_x)
%         ylabel('track mean (x)')
% 
%         subplot(3,2,3)
%         plot(times, trace_IKKa)
%         xlabel('Time')
%         ylabel('Intermediate (IKKa)')
%         ylim([0, max(trace_IKKa)])
% 
%         subplot(3,2,4)
%         plot(times, trace_I)
%         xlabel('Time')
%         ylabel('I')
%         
%         subplot(3,2,5)
%         plot(times, trace_Im)
%         xlabel('Time')
%         ylabel('Im')
%         
%         subplot(3,2,6)
%         plot(times, trace_NFkB)
%         xlabel('Time')
%         ylabel('NFkB')
%         ylim([0, 1])
%     end
    
%     if flag_plot == true
%         adap_out_color = 'magenta';
%         
%         figure()
%         times_in_s = times*60;
%         plot(times_in_s, trace_y, 'Color', adap_out_color, 'LineWidth', 2)
%         xlabel('time (s)')
%         ylabel('Output')
%         set(gca, 'fontsize', 18)
%     end
    
    if flag_plot == true
        input_color = 'blue';
        adap_out_color = 'magenta';
        
        figure()
        xlabel('time (min)')
        set(gca, 'fontsize', 18)
        
        yyaxis left
        plot(times, inputs, 'Color', input_color, 'LineWidth', 2)
        set(gca,'YColor',input_color)
        ylabel('Input')
        
        yyaxis right
        plot(times, trace_y, 'Color', adap_out_color, 'LineWidth', 2)
        set(gca,'YColor',adap_out_color)
        ylabel('Output')
        ylim([0, 1.05*max(trace_y)])     
        xlim([final_time/2, final_time/2+3*period])
    end
    
    if flag_plot == true
        figure('rend','painters','pos',[0 0 500 600])
        input_color = 'blue';
        adap_out_color = 'magenta';
        
        subplot(2, 1, 1)
        set(gca, 'fontsize', 22)
        plot(times, inputs, 'Color', input_color, 'LineWidth', 4)
        xlabel('time (\tau_0)')
        ylabel('Input')
        xlim([0, 2*period])
        set(findall(gcf,'-property','FontSize'),'FontSize',20)
        set(gca,'box','off')
        
        subplot(2, 1, 2)
        xlabel('time (\tau_0)')
        set(gca, 'fontsize', 22)
        plot(times, trace_y, 'Color', adap_out_color, 'LineWidth', 4)
        xlabel('time (\tau_0)')
        ylabel('Y')
        ylim([0, 1.05*max(trace_y)])     
        xlim([0, 2*period])
        set(findall(gcf,'-property','FontSize'),'FontSize',20)
        set(gca,'box','off')
    end

    if flag_plot == true
        my_ylim = 0.5;
        input_color = 'blue';
        NF_kB_color = 'red';
        figure()
        set(gca, 'fontsize', 18)
        if noise_str > 0
            title('Activation due to noise')
        else
            title('No activation without noise')
        end
        
        xlabel('time (min)')
        yyaxis left
        plot(times, inputs, 'LineWidth', 2, 'Color', input_color)
        ylim([0, my_ylim])
        set(gca,'YColor',input_color)
        ylabel('Input')
        yyaxis right
        plot(times, trace_NFkB, '-', 'LineWidth', 2, 'Color', NF_kB_color)
        ylim([0, 1])
        ylabel('NF-kB')
        set(gca,'YColor',NF_kB_color)
        set(gca,'box','off')
    end
    
    for i = 1:ensemblesize        
        activate_idx = find(trace_NFkB(:,i) > activation_threshold, 1);
        if isempty(activate_idx)
            time_activate(i) = times(end);
        else
            time_activate(i) = times(activate_idx);
        end
    end
    
    time_activate_avg = mean(time_activate);
    frac_activate = sum(time_activate < final_time)/ensemblesize;
    
    [~, peak_idx] = max(trace_NFkB, [], 1);
    time_peak = (peak_idx-1)*dt;
    
    
%     
%     index_half = floor(length(times)/2);
%     times   = times(index_half:end);
%     trace_y = trace_y(index_half:end);
%     Y_avg = trapz(times, trace_y)/(times(end) - times(1));

end
