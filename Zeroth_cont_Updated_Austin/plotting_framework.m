function plotting_framework(num_mon_A, num_mon_B, num_mon_C, num_input_A, num_input_B, num_input_C, u_sim_plot, q_svc_hat_plot,  PVGen_plot, loadprofile_plot, gamma_plot, f_obj_plot, nt, Plot_freq)
%initplot Plot to simplify the plotting framework

    num_mon = num_mon_A + num_mon_B + num_mon_C;
    num_input = num_input_A + num_input_B + num_input_C;

    % Initialize the line colors.
    line_colors = jet(max([num_mon, num_input]));

    % Data extraction
    u_sim_plotA = u_sim_plot(1:num_mon_A, :);
    u_sim_plotB = u_sim_plot((num_mon_A+1):(num_mon_A+num_mon_B), :);
    u_sim_plotC = u_sim_plot((num_mon_A+1):(num_mon_A+num_mon_B), :);

    q_svc_hat_plotA = q_svc_hat_plot(1:num_input_A, :);
    q_svc_hat_plotB =q_svc_hat_plot((num_input_A+1):(num_input_A+num_input_B), :);
    q_svc_hat_plotC =q_svc_hat_plot((num_input_A+num_input_B+1):num_input, :);

    % Set the figure size to 80% of the screen size
    %% Plotting
    %% PV and load
    subplot(5,2,1); 
    if nt > Plot_freq
        plot(nt-Plot_freq:nt, PVGen_plot, 'Color', 'r', 'LineWidth', 1); hold on;
    else
        plot(nt-Plot_freq+1:nt, PVGen_plot, 'Color', 'r', 'LineWidth', 1); hold on;
    end
    xlabel('Time (seconds)');
    ylabel('PV Gen Multiplier');
    grid on;

    subplot(5,2,2); 
    if nt > Plot_freq
        plot(nt-Plot_freq:nt, loadprofile_plot, 'Color', 'b', 'LineWidth', 1); hold on;
    else
        plot(nt-Plot_freq+1:nt, loadprofile_plot, 'Color', 'b', 'LineWidth', 1); hold on;
    end
    xlabel('Time (seconds)');
    ylabel('Load Multiplier');
    grid on;


    %% voltage with control
    subplot(5,2,3);
    if nt > Plot_freq
        for i = 1:num_mon_A
            % Plot the line.
            plot(nt-Plot_freq:nt, u_sim_plotA(i, :), 'Color', line_colors(i, :), 'LineWidth', 1); hold on;
        end
    else
        for i = 1:num_mon_A
            % Plot the line.
            plot(nt-Plot_freq+1:nt, u_sim_plotA(i, :), 'Color', line_colors(i, :), 'LineWidth', 1); hold on;
        end
    end
    ylim([0.9 1.1]);
    yticks([0.90 0.94 0.98 1.02 1.06 1.10]);
    xlabel('Time (seconds)');
    ylabel('Phase Voltage');
    title('Phase A with Control');
    grid on;
    
    subplot(5,2,5); 
    if nt > Plot_freq
        for i = 1:num_mon_B
            % Plot the line.
            plot(nt-Plot_freq:nt, u_sim_plotB(i, :), 'Color', line_colors(i, :), 'LineWidth', 1); hold on;
        end
    else
        for i = 1:num_mon_B
            % Plot the line.
            plot(nt-Plot_freq+1:nt, u_sim_plotB(i, :), 'Color', line_colors(i, :), 'LineWidth', 1); hold on;
        end
    end
    ylim([0.9 1.1]);
    yticks([0.90 0.94 0.98 1.02 1.06 1.10]);
    xlabel('Time (seconds)');
    ylabel('Phase Voltage');
    title('Phase B with Control');
    grid on;
    
    subplot(5,2,7); 
    if nt > Plot_freq
        for i = 1:num_mon_C
            % Plot the line.
            plot(nt-Plot_freq:nt, u_sim_plotC(i, :), 'Color', line_colors(i, :), 'LineWidth', 1); hold on;
        end
    else
        for i = 1:num_mon_C
            % Plot the line.
            plot(nt-Plot_freq+1:nt, u_sim_plotC(i, :), 'Color', line_colors(i, :), 'LineWidth', 1); hold on;
        end
    end
    ylim([0.9 1.1]);
    yticks([0.90 0.94 0.98 1.02 1.06 1.10]);
    xlabel('Time (seconds)');
    ylabel('Phase Voltage');
    title('Phase C with Control');
    grid on;
    
    %% reactive power with control
    subplot(5,2,4); 
    if nt > Plot_freq
        for i = 1:num_input_A
            % Plot the line.
            plot(nt-Plot_freq:nt, q_svc_hat_plotA(i, :), 'Color', line_colors(i, :), 'LineWidth', 1); hold on;
        end
    else
        for i = 1:num_input_A
            % Plot the line.
            plot(nt-Plot_freq+1:nt, q_svc_hat_plotA(i, :), 'Color', line_colors(i, :), 'LineWidth', 1); hold on;
        end
    end
    xlabel('Time (seconds)');
    ylabel('Q');
    title('Phase A Control');
    grid on;
    
    subplot(5,2,6); 
    if nt > Plot_freq
        for i = 1:num_input_B
            % Plot the line.
            plot(nt-Plot_freq:nt, q_svc_hat_plotB(i, :), 'Color', line_colors(i, :), 'LineWidth', 1); hold on;
        end
    else
        for i = 1:num_input_B
            % Plot the line.
            plot(nt-Plot_freq+1:nt, q_svc_hat_plotB(i, :), 'Color', line_colors(i, :), 'LineWidth', 1); hold on;
        end
    end
    xlabel('Time (seconds)');
    ylabel('Q');
    title('Phase B Control');
    grid on;
    
    subplot(5,2,8); 
    if nt > Plot_freq
        for i = 1:num_input_C
            % Plot the line.
            plot(nt-Plot_freq:nt, q_svc_hat_plotC(i, :), 'Color', line_colors(i, :), 'LineWidth', 1); hold on;
        end
    else
        for i = 1:num_input_C
            % Plot the line.
            plot(nt-Plot_freq+1:nt, q_svc_hat_plotC(i, :), 'Color', line_colors(i, :), 'LineWidth', 1); hold on;
        end
    end
    xlabel('Time (seconds)');
    ylabel('Q');
    title('Phase C Control');
    grid on;


    subplot(5,2,9); 
    if nt > Plot_freq
        for i = 1:num_input
            % Plot the line.
            semilogy(nt-Plot_freq:nt, (abs(gamma_plot(i, :))) , 'Color', line_colors(i, :), 'LineWidth', 1); hold on;
        end
    else
        for i = 1:num_input
            % Plot the line.
            semilogy(nt-Plot_freq+1:nt, (abs(gamma_plot(i, :))), 'Color', line_colors(i, :), 'LineWidth', 1); hold on;
        end
    end
    xlabel('Time (seconds)');
    ylabel('|\nabla f|');
    title('Gradient of the objective');
    grid on;

    subplot(5,2,10); 
    if nt > Plot_freq
        % Plot the line.
        semilogy(nt-Plot_freq:nt, abs(f_obj_plot), 'Color', 'r', 'LineWidth', 1); hold on;
    else
        % Plot the line.
        semilogy(nt-Plot_freq+2:nt, abs(f_obj_plot(2:nt)), 'Color', 'r', 'LineWidth', 1); hold on;
    end
    xlabel('Time (seconds)');
    ylabel('|f(t)|');
    title('Objective function');
    grid on;

end