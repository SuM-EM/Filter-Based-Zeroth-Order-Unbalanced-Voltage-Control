function [xi_new, mu_new, la_1_new, la_2_new, f_obj, gamma_new, grad_new, q_svc_new, q_svc_hat_new] = ...
    Cont_Zero_Dist(control, nt, u_sim_old, mu_old, xi_old, la_1_old, la_2_old, q_svc_max, q_svc_min, f_obj_old, gamma_old, q_svc_old, q_svc_old_2, num_mon, num_input)
%Controller The controller goes here. Controller do not control the main
%simulator.
%   Input:  ActNumBuses: Actual number of buses
%           ContCap: Controller rating (kW/KVar)
%           q_hat_old: Reactive power from controller in the previous step
%           v_hat: System-wide voltage at the previous time step
%           tt:  Current system time (for initial set point)
%   Output: q_hat: Desired controller injection

%% Run Simulation 
%% ------ Controller Parameter Setting ------------------------------%%
q_scale = 1e3;      u_scale = 1e3;
v_nom = 1e3;
u_nom = ones(num_mon,1)*v_nom;

%% Assumption 4 prime picking
x = primes(1000); % pick primes smaller than, say 1000
x = x(x>4);       % pick the ones larger than 4
delta = 0.1*(2+5*x(1:num_input)/x(num_input));

%% Controller 2
if control == 1
    c_q = 0.01; % quadratic control cost
    kappa = delta';
    kappa = kappa*0.05/0.5;

    amp = 0.01; % amplitude of perturbation signal = 0.2, 0.3 are best, 0.4 larger ossication, 0.1 is even larger
    dt = 2;    % unit:second, as discussed in main file

    eta = 1e-7;  % 1e-4 is too large causing ossication
end

%% Controller 3
if control == 2
    c_q = 0.01; % quadratic control cost
    kappa = delta';
    kappa = kappa*0.05/0.5;

    amp = 0.01;    % amplitude of perturbation signal = 0.2, 0.3 are best, 0.4 larger ossication, 0.1 is even larger
    dt = 2;       % unit:second, as discussed in main file

    eta    = 0.1*1e-7;  %0.0001 was
    alpha  = 0.6;
end

if control == 1
    %% Update Control  -- Controller 2 - Projected Gradient =====================================================
    if nt > 2
        f_obj = sum((u_scale*u_sim_old.^2 -u_nom).^2);
        gamma_new = 2/amp*f_obj*sin(2*pi*kappa*(nt*dt)) + c_q*2*q_svc_old;
        q_svc_new = max(min(q_svc_old - eta*gamma_new, q_svc_max/q_scale),q_svc_min/q_scale);
    elseif nt == 2
        f_obj = sum((u_scale*u_sim_old.^2 -u_nom).^2);
        gamma_new = gamma_old;
        q_svc_new = q_svc_old;
    else
        % Initialize all variable to zero
        f_obj = f_obj_old;
        gamma_new = gamma_old;
        q_svc_new = q_svc_old;
    end
    xi_new = xi_old;        mu_new = mu_old;
    la_1_new = la_1_old;    la_2_new = la_2_old;      grad_new = gamma_new;
elseif control == 2
    %% Update Control  -- Controller 2 - Projected Gradient with lpf-hpf ========================================
    if nt > 2
        f_obj = sum((u_scale*u_sim_old.^2 -u_nom).^2);
        gamma_new = (f_obj - f_obj_old);
        grad = 2/amp*gamma_new.*sin(2*pi*kappa*(nt*dt)) - alpha/eta*(q_svc_old-q_svc_old_2);
        q_svc_new = max(min(q_svc_old - eta*(2/amp*gamma_new.*sin(2*pi*kappa*(nt*dt)) + 2*c_q*q_svc_old) + alpha*(q_svc_old-q_svc_old_2)...
                                                                                            , q_svc_max/q_scale), q_svc_min/q_scale);
    elseif nt == 2
        f_obj = sum((u_scale*u_sim_old.^2 -u_nom).^2) + c_q*sum(q_svc_old.^2);
        gamma_new = gamma_old;
        grad = 2/amp*gamma_new.*sin(2*pi*kappa*(nt*dt));
        q_svc_new = q_svc_old;
    else
        f_obj = f_obj_old;
        gamma_new = gamma_old;
        grad = 2/amp*gamma_new.*sin(2*pi*kappa*(nt*dt));
        q_svc_new = q_svc_old;
    end
    xi_new = xi_old;        mu_new = mu_old;
    la_1_new = la_1_old;    la_2_new = la_2_old;      grad_new = grad;
end

    % 1- Add Perturbation
    q_svc_hat_new = q_scale*(q_svc_new + amp*sin(2*pi*kappa*((nt-1)*dt)));

end

