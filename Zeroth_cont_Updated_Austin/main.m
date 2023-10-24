clc; clear all; close all; beep off;

%% Missing string for master.dss
soucepath  = pwd;
dssfile    = '\syn-Austin-TDgrid-v03\syn-austin-D_only-v03\P1U\base\opendss\p1uhs17_1247\p1uhs17_1247--p1udt1947\Master.dss';
% dssfile    = '\13Bus\IEEE13Nodeckt.dss';

%% Identify active buses
[ArrangedBusName,phasesinbuses,ActNumBuses] = opendss_startup(dssfile,soucepath);

% PV Nodes == p1udt22135.1.2.3 p1udt1947.1.2.3 p1udt24439.1.2.3  p1udt24957.1.2.3  p1udt3662.1.2.3    p1udt23231.1.2.3
%   p1udt3660.1.2.3  p1udt19457.1.2.3  p1udt19634.1.2.3   p1udt24549.1.2.3
%   p1udt18509.1.2.3 p1udt21176.1.2.3   p1udt19633.1.2.3  p1udt10925.1.2.3
%    p1udt13860.1.2.3  p1udt20754.1.2.3  p1udt22990.1.2.3
%    p1udt24247.1.2.3  p1udt24341.1.2.3  p1umv27.1.2.3

%% Select PV locations; locations are selected based on availability of phase in a particular bus
PV_cap = 200;
PVLocA = {'p1udt22135.1', 'p1udt1947.1', 'p1udt24439.1', 'p1udt24957.1', 'p1udt3662.1', 'p1udt23231.1', 'p1udt3660.1', 'p1udt19457.1', 'p1udt19634.1', 'p1udt24549.1', 'p1udt18509.1', 'p1udt21176.1', 'p1udt19633.1', 'p1udt10925.1', 'p1udt13860.1', 'p1udt20754.1', 'p1udt22990.1', 'p1udt24247.1', 'p1udt24341.1', 'p1umv27.1'};
PVLocB = {'p1udt22135.2', 'p1udt1947.2', 'p1udt24439.2', 'p1udt24957.2', 'p1udt3662.2', 'p1udt23231.2', 'p1udt3660.2', 'p1udt19457.2', 'p1udt19634.2', 'p1udt24549.2', 'p1udt18509.2', 'p1udt21176.2', 'p1udt19633.2', 'p1udt10925.2', 'p1udt13860.2', 'p1udt20754.2', 'p1udt22990.2', 'p1udt24247.2', 'p1udt24341.2', 'p1umv27.2'};
PVLocC = {'p1udt22135.3', 'p1udt1947.3', 'p1udt24439.3', 'p1udt24957.3', 'p1udt3662.3', 'p1udt23231.3', 'p1udt3660.3', 'p1udt19457.3', 'p1udt19634.3', 'p1udt24549.3', 'p1udt18509.3', 'p1udt21176.3', 'p1udt19633.3', 'p1udt10925.3', 'p1udt13860.3', 'p1udt20754.3', 'p1udt22990.3', 'p1udt24247.3', 'p1udt24341.3', 'p1umv27.3'};
%% in Phase A
location = zeros(ActNumBuses,1); % Initialize PV locations
for ii = 1:numel(PVLocA)
    location(strcmp(ArrangedBusName(:,1), PVLocA{ii})) = 1;
end
capacity = PV_cap*ones(ActNumBuses,1);     % Select PV capacity array
PVinPhaseA = phasesinbuses(:,1).*location;   
PVCapPhaseA = PVinPhaseA.*capacity;
%% in Phase B
location = zeros(ActNumBuses,1); % Initialize PV locations
for ii = 1:numel(PVLocB)
    location(strcmp(ArrangedBusName(:,2), PVLocB{ii})) = 1;
end
capacity = PV_cap*ones(ActNumBuses,1);     % Select PV capacity array
PVinPhaseB = phasesinbuses(:,2).*location;  
PVCapPhaseB = PVinPhaseB.*capacity;
%% in Phase C
location = zeros(ActNumBuses,1); % Initialize PV locations
for ii = 1:numel(PVLocC)
    location(strcmp(ArrangedBusName(:,3), PVLocC{ii})) = 1;
end
capacity = PV_cap*ones(ActNumBuses,1);     % Select PV capacity array
PVinPhaseC = phasesinbuses(:,3).*location; 
PVCapPhaseC = PVinPhaseC.*capacity;

%% Capture PV and Load Profile
% T = readtable([soucepath,'\15minute_data_austin\15minute_data_austin.csv'],'Range','A1:CA1000');  
% loadprofile = (table2array(T(708:803,{'grid'}))+table2array(T(708:803,{'solar'})))/7;
% PVData = table2array(T(708:803,{'solar'}))/5;
% loadunimul = 4.0; %% Load multiplier for system-wide load increment

T = readtable([soucepath,'\15minute_data\15minute_data.csv'],'Range','A1:CA100');  
% Note that the multipliers are now removed
loadprofile = (2/4)*(table2array(T(1:96,{'grid'}))+table2array(T(1:96,{'solar'})));
PVData = abs(table2array(T(1:96,{'solar'})))/5; PVData(PVData < 1e-2) = 0;
loadunimul = 1.00;                      % Load multiplier for system-wide load increment
% loadunimul = 4.0;                      % Load multiplier (old, load profile multiplier 1/7, pv data multiplier 1/5)
% loadprofile = ones(1,96);
% PVData = zeros(1,96);
% loadunimul = 1.0; %% Load multiplier for system-wide load increment
T_hour = 0:0.25:(24-0.25);

%% Select Controller locations; locations are selected based on availability of phase in a particular bus
% Monitoring Nodes == p1udt19983.1.2.3    p1udt22005.1.2.3  p1udt15201.1.2.3  p1udt7267.1.2.3   p1udt16652.1.2.3   p1udt24957.1.2.3  p1udt110.1.2.3  p1udt15202.1.2.3
% Control Nodes == p1udt24141.1.2.3  p1udt22383.1.2.3   p1udt13282.1.2.3  p1udt23806.1.2.3  p1udt15201.1.2.3  p1udt3662.1.2.3   p1udt19457.1.2.3  p1udt19633.1.2.3  p1udt7267.1.2.3
%% Same bus for monitoring and control action
SVC_cap = 10000; % in KVAr
MonLocA = {'p1udt19983.1', 'p1udt22005.1', 'p1udt15201.1', 'p1udt7267.1', 'p1udt16652.1', 'p1udt24957.1', 'p1udt110.1', 'p1udt15202.1'};
MonLocB = {'p1udt19983.2', 'p1udt22005.2', 'p1udt15201.2', 'p1udt7267.2', 'p1udt16652.2', 'p1udt24957.2', 'p1udt110.2', 'p1udt15202.2'};
MonLocC = {'p1udt19983.3', 'p1udt22005.3', 'p1udt15201.3', 'p1udt7267.3', 'p1udt16652.3', 'p1udt24957.3', 'p1udt110.3', 'p1udt15202.3'};
ControlLocA = {'p1udt24141.1', 'p1udt22383.1', 'p1udt13282.1', 'p1udt23806.1', 'p1udt15201.1', 'p1udt3662.1', 'p1udt19457.1', 'p1udt19633.1', 'p1udt7267.1'};
ControlLocB = {'p1udt24141.2', 'p1udt22383.2', 'p1udt13282.2', 'p1udt23806.2', 'p1udt15201.2', 'p1udt3662.2', 'p1udt19457.2', 'p1udt19633.2', 'p1udt7267.2'};
ControlLocC = {'p1udt24141.3', 'p1udt22383.3', 'p1udt13282.3', 'p1udt23806.3', 'p1udt15201.3', 'p1udt3662.3', 'p1udt19457.3', 'p1udt19633.3', 'p1udt7267.3'};

%% Monitoring
MonindexA = []; MonindexB = []; MonindexC = [];
for ii = 1:numel(MonLocA)
    MonindexA(ii) = find(strcmp(ArrangedBusName(:,1), MonLocA{ii}));
end
for ii = 1:numel(MonLocB)
    MonindexB(ii) = find(strcmp(ArrangedBusName(:,2), MonLocB{ii}));
end
for ii = 1:numel(MonLocC)
    MonindexC(ii) = find(strcmp(ArrangedBusName(:,3), MonLocC{ii}));
end

%% Control
indexA = []; indexB = []; indexC = [];
%% in Phase A
location = zeros(ActNumBuses,1); % Randomly select PV locations
for ii = 1:numel(ControlLocA)
    indexA(ii) = find(strcmp(ArrangedBusName(:,1), ControlLocA{ii}));
    location(indexA(ii)) = 1;
end
capacity = SVC_cap*ones(ActNumBuses,1);     % Randomly select PV capacity
ContinPhaseA = phasesinbuses(:,1).*location';  
ContCapPhaseA = ContinPhaseA.*capacity;
%% in Phase B
location = zeros(ActNumBuses,1); % Randomly select PV locations
for ii = 1:numel(ControlLocB)
    indexB(ii) = find(strcmp(ArrangedBusName(:,2), ControlLocB{ii}));
    location(indexB(ii)) = 1;
end
capacity = SVC_cap*ones(ActNumBuses,1);     % Randomly select PV capacity
ContinPhaseB = phasesinbuses(:,2).*location';  
ContCapPhaseB = ContinPhaseB.*capacity;
%% in Phase C
location = zeros(ActNumBuses,1); % Randomly select PV locations
for ii = 1:numel(ControlLocC)
    indexC(ii) = find(strcmp(ArrangedBusName(:,3), ControlLocC{ii}));
    location(indexC(ii)) = 1;
end
capacity = SVC_cap*ones(ActNumBuses,1);     % Randomly select PV capacity
ContinPhaseC = phasesinbuses(:,3).*location';  
ContCapPhaseC = ContinPhaseC.*capacity;

ContCap = [ContCapPhaseA, ContCapPhaseB, ContCapPhaseC]; 

%% Tap Positions
% OLTC
TAPOLTC = [1; 1; 1];
%SVR
TAPSVR = [1; 1; 1];

%% Controller initialization
% Interpolation of data for faster dynamics
T_sec = 0:(1/(5*360)):(24-0.25); %1/5 sec of gap between two time steps
loadprofile_sec_2 = interp1(T_hour,loadprofile, T_sec,'spline');
loadprofile_sec = movmean(loadprofile_sec_2,200);
PVData_sec_2 = interp1(T_hour,PVData, T_sec,'spline'); PVData_sec_2(PVData_sec_2<1e-3) = 0;
PVData_sec = movmean(PVData_sec_2,200);

% select time period for simulation
T_index_sim = 1:5*360*(24-0.25);       
T_sim = T_sec(T_index_sim);

PVData = PVData_sec(T_index_sim);
loadprofile = loadprofile_sec(T_index_sim);

num_sim_time = length(T_sim);

% Monitored Nodes
bus_volt_A = MonindexA;
bus_volt_B = MonindexB; 
bus_volt_C = MonindexC; 
num_mon_A = length(bus_volt_A);
num_mon_B = length(bus_volt_B);
num_mon_C = length(bus_volt_C);
num_mon = num_mon_A + num_mon_B + num_mon_C;

% Zeroth-Order Controller nodes
bus_ctrl_A = indexA;
bus_ctrl_B = indexB;
bus_ctrl_C = indexC;
num_input_A = length(bus_ctrl_A);
num_input_B = length(bus_ctrl_B);
num_input_C = length(bus_ctrl_C);
num_input = num_input_A + num_input_B + num_input_C;

% Controller parameter
SVC_cap = SVC_cap;  %Unit: MVAR
q_svc_min = - SVC_cap*ones(num_input,1); 
q_svc_max =  SVC_cap*ones(num_input,1); 

% Controller Initialization
q_hat = zeros(ActNumBuses,3,numel(PVData));   % Reactive power injected
v_hat = zeros(ActNumBuses,3,numel(PVData));   % system wide voltages calculated

% voltage variables
u_sim = zeros(num_mon,num_sim_time);

% control variables
q_svc = zeros(num_input, num_sim_time);
q_svc_hat = zeros(num_input, num_sim_time);

% dual variables
la_1 = zeros(num_mon, num_sim_time);
la_2 = zeros(num_mon, num_sim_time);

% fast variables
xi = zeros(num_mon,num_input,num_sim_time);
mu = zeros(num_mon,num_sim_time);
mu(:,1) = ones(num_mon,1);

% Controller 2
f_obj = zeros(1, num_sim_time);
gamma = zeros(num_input, num_sim_time);     
grad = zeros(num_input, num_sim_time);       % This is extracted gradient

%% Controller Select
control = 2;   %% Choose between controller 1 and 2

%% Initial plotting structure
% Let us plot only two of the monitoed node and two of the control node data
num_mon_A_plot = 2; num_mon_B_plot = 2; num_mon_C_plot = 2;
num_input_A_plot = 2; num_input_B_plot = 2; num_input_C_plot = 2;
num_mon_plot = num_mon_A_plot + num_mon_B_plot + num_mon_C_plot;
num_input_plot = num_input_A_plot + num_input_B_plot + num_input_C_plot;

u_sim_plot = zeros(num_mon_plot,2);
q_svc_hat_plot = zeros(num_input_plot,2);
grad_plot = zeros(num_input_plot, 2);
f_obj_plot = zeros(1, 2);
PVGen_plot = [0 0];
loadprofile_plot = [0 0];
Plot_freq = 50;   % How frequently do you want the plot to refresh?
figure('Units', 'Normalized', 'Position', [0 0 0.8 0.8]);
plotting_framework(num_mon_A_plot, num_mon_B_plot, num_mon_C_plot, num_input_A_plot, num_input_B_plot, num_input_C_plot, u_sim_plot, q_svc_hat_plot, PVGen_plot, loadprofile_plot, grad_plot, f_obj_plot, 2, 2); 
ntinc = 1;  % Initialize plotting pointer

%% Run Simulation 
for nt = 1:numel(PVData)
    
    if nt > 2
        u_sim_old = u_sim(:,nt-1);
	    mu_old = mu(:,nt-1);            xi_old = xi(:,:,nt-1);
	    la_1_old = la_1(:,nt-1);        la_2_old = la_2(:,nt-1); 
        q_svc_old = q_svc(:,nt-1);      q_svc_old_2 = q_svc(:,nt-2);

        f_obj_old = f_obj(:,nt-1);
        gamma_old = gamma(:,nt-1);    
    elseif nt == 2
        u_sim_old = u_sim(:,nt-1);
	    mu_old = mu(:,nt-1);            xi_old = xi(:,:,nt-1);
	    la_1_old = la_1(:,nt-1);        la_2_old = la_2(:,nt-1); 
        q_svc_old = q_svc(:,nt-1);      q_svc_old_2 = q_svc(:,nt-1);

        f_obj_old = f_obj(:,nt-1);
        gamma_old = gamma(:,nt-1);

    else
        % -- Initialization is needed
        u_sim_old = u_sim(:,1);
	    mu_old = mu(:,1);            xi_old = xi(:,:,1);
	    la_1_old = la_1(:,1);        la_2_old = la_2(:,1); 
        q_svc_old = q_svc(:,1);      q_svc_old_2 = q_svc(:,1);

        f_obj_old = f_obj(:,1);
        gamma_old = gamma(:,1);
    end

    %% Call the controller (note time control is outside)
    [xi_new, mu_new, la_1_new, la_2_new, f_obj_val, gamma_new, grad_new, q_svc_new, q_svc_hat_new] = ...
                Cont_Zero_Dist(control, nt, u_sim_old, mu_old, xi_old, la_1_old, la_2_old, q_svc_max, q_svc_min, f_obj_old, gamma_old, q_svc_old, q_svc_old_2, num_mon, num_input);

    xi(:,:,nt) = xi_new;            mu(:,nt) = mu_new; 	
	la_1(:,nt) = la_1_new;          la_2(:,nt) = la_2_new; 
    q_svc(:,nt) = q_svc_new;        q_svc_hat(:,nt) = q_svc_hat_new;

    f_obj(:,nt) = f_obj_val;
    gamma(:,nt) = gamma_new;        grad(:,nt) = grad_new;

    % 2- reformat the controller
    q_hat(bus_ctrl_A,1,nt) = q_svc_hat(1:num_input_A,nt); % phase A 
    q_hat(bus_ctrl_B,2,nt) = q_svc_hat((num_input_A+1):(num_input_A+num_input_B),nt); % phase B
    q_hat(bus_ctrl_C,3,nt) = q_svc_hat((num_input_A+num_input_B+1):num_input,nt); % phase C


    %% Pre-process data for simulator
    PVGen = PVData(nt);
    LoadMul = loadunimul*loadprofile(nt);
    PVCap = [PVCapPhaseA,PVCapPhaseB,PVCapPhaseC];
    QINJ = round(q_hat(:,:,nt)*1000)/1000;

    %% Processed TAP Positions
    TAPPOS = [TAPOLTC; TAPSVR];
    if nt >= 6900 && nt <=9000
        TAPPOS = 2*[TAPOLTC; TAPSVR];
    elseif nt >=14000 && nt <17500
        TAPPOS = 0*[TAPOLTC; TAPSVR]; TAPPOS(2:3) = -3*TAPOLTC(2:3); TAPPOS(5:6) = -3*TAPSVR(2:3);
    elseif nt >=17500 && nt <20000
        TAPPOS = -1*[TAPOLTC; TAPSVR]; TAPPOS(2:3) = -6*TAPOLTC(2:3); TAPPOS(5:6) = -6*TAPSVR(2:3);
    elseif nt >=20000 && nt <25000
        TAPPOS = -1*[TAPOLTC; TAPSVR]; TAPPOS(2:3) = -8*TAPOLTC(2:3); TAPPOS(5:6) = -8*TAPSVR(2:3);
    elseif nt >=25000 && nt <27500
        TAPPOS = -1*[TAPOLTC; TAPSVR]; TAPPOS(2:3) = -6*TAPOLTC(2:3); TAPPOS(5:6) = -6*TAPSVR(2:3);
    elseif nt >=27500 && nt <30000
        TAPPOS = 0*[TAPOLTC; TAPSVR]; TAPPOS(2:3) = -3*TAPOLTC(2:3); TAPPOS(5:6) = -3*TAPSVR(2:3);
    end
%     TAPPOS = 0*[TAPOLTC; TAPSVR];   % Control if no SVR are there

    %% Call Simulator and store voltage
    [ArrangedVoltagesA,ArrangedVoltagesB,ArrangedVoltagesC] = ...
                    opendss_interface(PVGen,LoadMul,ArrangedBusName,PVCap,ContCap,QINJ,TAPPOS,phasesinbuses,ActNumBuses,dssfile,soucepath);
    
    v_hat(:,:,nt) = [ArrangedVoltagesA,ArrangedVoltagesB,ArrangedVoltagesC];

    u_sim(1:num_mon_A, nt) = ArrangedVoltagesA(bus_volt_A);
    u_sim((num_mon_A+1):(num_mon_A+num_mon_B), nt) = ArrangedVoltagesB(bus_volt_B);
    u_sim((num_mon_A+num_mon_B+1):num_mon, nt) = ArrangedVoltagesC(bus_volt_C);
    
    %% Process voltage for live nodes
    ArrangedVoltagesA(phasesinbuses(:,1)==0) = []; 
    ArrangedVoltagesPhA_Time(:,nt) = ArrangedVoltagesA;
    ArrangedVoltagesB(phasesinbuses(:,2)==0) = [];
    ArrangedVoltagesPhB_Time(:,nt) = ArrangedVoltagesB;
    ArrangedVoltagesC(phasesinbuses(:,3)==0) = [];
    ArrangedVoltagesPhC_Time(:,nt) = ArrangedVoltagesC;

    Qinjected(:,:,nt) = QINJ;

    %% Plotting related variables  --> Deleted latter
    q_svc_hat_plot(:,ntinc) = q_svc_hat([1:num_input_A_plot,(num_input_A+1):(num_input_A+num_input_B_plot),(num_input_A+num_input_B+1):(num_input_A+num_input_B+num_input_C_plot)],nt);
    u_sim_plot(:,ntinc) = u_sim([1:num_mon_A_plot,(num_mon_A+1):(num_mon_A+num_mon_B_plot),(num_mon_A+num_mon_B+1):(num_mon_A+num_mon_B+num_mon_C_plot)], nt);
    PVGen_plot(ntinc) = PVData(nt);
    loadprofile_plot(ntinc) = loadunimul*loadprofile(nt);
    grad_plot(:,ntinc) = grad([1:num_input_A_plot,(num_input_A+1):(num_input_A+num_input_B_plot),(num_input_A+num_input_B+1):(num_input_A+num_input_B+num_input_C_plot)], nt);

    f_obj_plot(ntinc) = f_obj(nt);

    ntinc = ntinc + 1;  % plotting pointer update
   
    if rem(nt, Plot_freq) == 0      % Update plot after every 10 iteration
        if nt > Plot_freq
            q_svc_hat_plot = [q_svc_hat([1:num_input_A_plot,(num_input_A+1):(num_input_A+num_input_B_plot),(num_input_A+num_input_B+1):(num_input_A+num_input_B+num_input_C_plot)],nt-Plot_freq), q_svc_hat_plot];
            u_sim_plot = [u_sim([1:num_mon_A_plot,(num_mon_A+1):(num_mon_A+num_mon_B_plot),(num_mon_A+num_mon_B+1):(num_mon_A+num_mon_B+num_mon_C_plot)],nt-Plot_freq), u_sim_plot];
            PVGen_plot = [PVData(nt-Plot_freq), PVGen_plot];
            loadprofile_plot = [loadunimul*loadprofile(nt-Plot_freq), loadprofile_plot];
            grad_plot = [grad([1:num_input_A_plot,(num_input_A+1):(num_input_A+num_input_B_plot),(num_input_A+num_input_B+1):(num_input_A+num_input_B+num_input_C_plot)], nt-Plot_freq), grad_plot];
            f_obj_plot = [f_obj(nt-Plot_freq), f_obj_plot];
        end
        tic;
        pause(0.01)
        plotting_framework(num_mon_A_plot, num_mon_B_plot, num_mon_C_plot, num_input_A_plot, num_input_B_plot, num_input_C_plot, u_sim_plot, q_svc_hat_plot, PVGen_plot, loadprofile_plot, grad_plot, f_obj_plot, nt, Plot_freq)
        ntinc = 1;  % Initialize plotting pointer
        q_svc_hat_plot = []; u_sim_plot = []; PVGen_plot = []; loadprofile_plot =[]; grad_plot = []; f_obj_plot = [];
        toc
    end

end

writematrix(f_obj,'f_obj_control.txt')
writematrix(u_sim,'u_sim_control.txt')
writematrix(q_svc_hat,'q_svc_hat_control.txt')