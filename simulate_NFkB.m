function Y_avg = simulate_NFkB(offset, amplitude, period, noise_str, noise_corr_time, flag_plot, final_time)
    % Author: Weerapat Pittayakanchit
    % Simulate the incoherent feedforward from James E. Ferrel, Jr.'s paper
    % Perfect and Near-Perfect Adaptation
    
    %% Set the simulation specfic details
%     final_time = 500;
    t_initial = 20;
    dt = 0.0001;
    times = 0:dt:final_time;
    
    square_waves = square(times*(2*pi)/period);
    inputs = offset*ones(length(times), 1);
    noise = 0;
    t_next = 0;
    index = 1;
    while index <= length(times)
        time  = times(index);
        if (noise_str ~= 0) && (time >= t_next)
            t_next = t_next + exprnd(noise_corr_time);
            noise  = noise_str*(2*rand() - 1);
        end
            
        if time > t_initial
            inputs(index) = max(offset + amplitude*(1 + square_waves(index))/2 + noise, 0);
        else
            inputs(index) = max(offset + noise, 0);
        end
        index = index + 1;
    end
    %% Set the initial conditions of the variables in the model
    ensemblesize = 1;
%     [0.0903    0.4941    3.5348    0.1254    0.6993];
%     Nn   = zeros(ensemblesize, 1);
%     Im   = zeros(ensemblesize, 1);
%     I    = zeros(ensemblesize, 1);
%     IKKa = zeros(ensemblesize, 1);
%     IKKi = zeros(ensemblesize, 1);
%     Nn   = 0.0903*ones(ensemblesize, 1);
%     Im   = 0.4941*ones(ensemblesize, 1);
%     I    = 3.5348*ones(ensemblesize, 1);
%     IKKa = 0.1254*ones(ensemblesize, 1);
%     IKKi = 0.6993*ones(ensemblesize, 1);
%     [0.0785, 0.3737, 3.6224, 0.0936, 0.5020]
    Nn   = 0.0785*ones(ensemblesize, 1);
    Im   = 0.3737*ones(ensemblesize, 1);
    I    = 3.6224*ones(ensemblesize, 1);
    IKKa = 0.0936*ones(ensemblesize, 1);
    IKKi = 0.5020*ones(ensemblesize, 1);
    
    %% Set observable variables to be plotted
    trace_NFkB = zeros(length(times), 1);
    trace_IKKa = zeros(length(times), 1);
    trace_IKKi = zeros(length(times), 1);
    trace_I    = zeros(length(times), 1);
    trace_Im   = zeros(length(times), 1);
    
    %% Simulate
    for index = 1:length(times)
        time  = times(index);
        input = inputs(index);
        TNF   = input;
        
        %% Parameters of the models
        kNin = 5.4;         klin = 0.018;   kt = 1.03;
        ktl  = 0.24;        KI  = 0.035;    KN = 0.029; 
        gamma_m  = 0.017;   alpha = 1.05;   Ntot = 1;
        ka   = 0.24;        ki = 0.18;      kp = 0.036;
        ka20 = 0.0018;      IKKt = 2;       a20 = 0.0026;
    
        dNn     = dt*(  kNin * (Ntot - Nn) * KI./(KI + I)    - klin * I .* Nn ./(Nn + KN));
        dIm     = dt*(  kt * Nn.^2                           - gamma_m * Im);
        dI      = dt*(  ktl * Im                             - alpha * IKKa .* (Ntot - Nn) .* I ./(KI + I));
        dIKKa   = dt*(  ka .* TNF .* (IKKt - IKKa - IKKi)    - ki .* IKKa);
        dIKKi   = dt*(  ki * IKKa                            - kp * IKKi * ka20./(ka20 + a20 * TNF));
        
        Nn      = Nn + dNn;
        Im      = Im + dIm;
        I       = I  + dI;
        IKKa    = IKKa + dIKKa;
        IKKi    = IKKi + dIKKi;
      
        % track X and Y
        trace_NFkB(index) = Nn(1);
        trace_IKKa(index) = IKKa(1);
        trace_IKKi(index) = IKKi(1);
        trace_I(index)    = I(1);
        trace_Im(index)   = Im(1);
    end
    
    %% Plot

    if flag_plot == true
        figure()
        set(gca, 'fontsize', 18)
        subplot(3,2,1)
        plot(times, inputs)
        ylabel('Input')

        subplot(3,2,2)
        plot(times, trace_IKKi)
        ylabel('Feedback (IKKi)')

        subplot(3,2,3)
        plot(times, trace_IKKa)
        xlabel('Time')
        ylabel('Output (IKKa)')

        subplot(3,2,4)
        plot(times, trace_I)
        xlabel('Time')
        ylabel('I')
        
        subplot(3,2,5)
        plot(times, trace_Im)
        xlabel('Time')
        ylabel('Im')
        
        subplot(3,2,6)
        plot(times, trace_NFkB)
        xlabel('Time')
        ylabel('NFkB')
    end

%     if flag_plot == true
%         figure()
%         subplot(3, 1 ,1)
%         plot(times, inputs)
%         ylabel('Input')
% 
%         subplot(3,1,2)
%         plot(times, trace_IKKi)
%         ylabel('Feedback (IKKi)')
% 
%         subplot(3,1,3)
%         plot(times, trace_IKKa)
%         xlabel('Time')
%         ylabel('Output (IKKa)')
%     end
    
    index_half = floor(length(times)/2);
    times   = times(index_half:end);
    trace_IKKa = trace_IKKa(index_half:end);
    Y_avg = trapz(times, trace_IKKa)/(times(end) - times(1));

end
