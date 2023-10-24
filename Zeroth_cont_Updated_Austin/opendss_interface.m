function [ArrangedVoltagesA,ArrangedVoltagesB,ArrangedVoltagesC] = ...
                            opendss_interface(PVGen,LoadMul,ArrangedBusName,PVCapacity,ContCap,QINJ,TAPPOS,phasesinbuses,ActNumBuses,dssfile,soucepath)
%opendss_interface Input:  PVGen:              Global PV Generation
%                          ArrangedBusNameOld: Which buses has which phases active
%                          PVCapacity:         What is PV Capacity in each phases, node-wise
%                          ContCap:            What is controller capacity in buses
%                          QINJ:               Controller outputs
%                          phasesinbuses:      Which node has which phases
%                          ActNumBuses:        Actual number of buses
%                          masterfile:         A string representing missing path of austiin grid
%   Calculates System-wide voltages using OpenDSS
%                  Output: ArrangedVoltagesA: PU Voltages in Phase A
%                          ArrangedVoltagesB: PU Voltages in Phase B
%                          ArrangedVoltagesC: PU Voltages in Phase C

PVCapacityinPhaseA = PVCapacity(:,1);
PVCapacityinPhaseB = PVCapacity(:,2);
PVCapacityinPhaseC = PVCapacity(:,3);

ContCapacityinPhaseA = ContCap(:,1);
ContCapacityinPhaseB = ContCap(:,2);
ContCapacityinPhaseC = ContCap(:,3);


%% Run at the beginning of the session
DSS = DSS_MATLAB.IDSS;

%% Used aliases
Text = DSS.Text;
Circuit = DSS.ActiveCircuit;
Solution = DSS.ActiveCircuit.Solution;
Bus = DSS.ActiveCircuit.ActiveBus;
Load = DSS.ActiveCircuit.Loads;
DSS.ClearAll();

Text.Command = ['Redirect [', soucepath, dssfile,']'];

Text.Command=['set mode =snap loadmult=' num2str(LoadMul)];   %% Control overall loads here

%% Define all the PV outputs
PV_num = 1;
for ii = 1: ActNumBuses
    if PVCapacityinPhaseA(ii) > 1
        Text.Command = ['new Generator.PV_A_', num2str(PV_num), ' bus=', ArrangedBusName{ii,1}, ' Model=1 kw=', num2str(PVCapacityinPhaseA(ii)*PVGen), ' kvar=0'];
        PV_num = PV_num + 1;
    end
end

PV_num = 1;
for ii = 1: ActNumBuses
    if PVCapacityinPhaseB(ii) > 1
        Text.Command = ['new Generator.PV_B_', num2str(PV_num), ' bus=', ArrangedBusName{ii,2}, ' Model=1 kw=', num2str(PVCapacityinPhaseB(ii)*PVGen), ' kvar=0'];
        PV_num = PV_num + 1;
    end
end

PV_num = 1;
for ii = 1: ActNumBuses
    if PVCapacityinPhaseC(ii) > 1
        Text.Command = ['new Generator.PV_C_', num2str(PV_num), ' bus=', ArrangedBusName{ii,3}, ' Model=1 kw=', num2str(PVCapacityinPhaseC(ii)*PVGen), ' kvar=0'];
        PV_num = PV_num + 1;
    end
end

%% Define all the Controller outputs

Cont_num = 1;
for ii = 1: ActNumBuses
    if QINJ(ii,1) ~= 0
        Text.Command = ['new Generator.Cont_A_', num2str(Cont_num), ' bus=', ArrangedBusName{ii,1}, ' Model=1 kw=0 kvar=', num2str(QINJ(ii,1))];
        Cont_num = Cont_num + 1;
    end
end

Cont_num = 1;
for ii = 1: ActNumBuses
    if QINJ(ii,2) ~= 0
        Text.Command = ['new Generator.Cont_B_', num2str(Cont_num), ' bus=', ArrangedBusName{ii,2}, ' Model=1 kw=0 kvar=', num2str(QINJ(ii,2))];
        Cont_num = Cont_num + 1;
    end
end

Cont_num = 1;
for ii = 1: ActNumBuses
    if QINJ(ii,3) ~= 0
        Text.Command = ['new Generator.Cont_C_', num2str(Cont_num), ' bus=', ArrangedBusName{ii,3}, ' Model=1 kw=0 kvar=', num2str(QINJ(ii,3))];
        Cont_num = Cont_num + 1;
    end
end

%% Set Regulators
% OLTC
Text.Command = ['new regcontrol.Reg1  transformer=Reg1 winding=2  vreg=122  band=2  ptratio=20 ctprim=700  R=3   X=9  TapNum=', num2str(TAPPOS(1)), ' maxtapchange=0'];
Text.Command = ['new regcontrol.Reg2  transformer=Reg2 winding=2  vreg=122  band=2  ptratio=20 ctprim=700  R=3   X=9  TapNum=', num2str(TAPPOS(2)), ' maxtapchange=0'];
Text.Command = ['new regcontrol.Reg3  transformer=Reg3 winding=2  vreg=122  band=2  ptratio=20 ctprim=700  R=3   X=9  TapNum=', num2str(TAPPOS(3)), ' maxtapchange=0'];

% SVR
Text.Command = ['new regcontrol.Reg4  transformer=Reg4 winding=2  vreg=122  band=2  ptratio=20 ctprim=700  R=3   X=9 TapNum=', num2str(TAPPOS(4)), ' maxtapchange=0'];
Text.Command = ['new regcontrol.Reg5  transformer=Reg5 winding=2  vreg=122  band=2  ptratio=20 ctprim=700  R=3   X=9 TapNum=', num2str(TAPPOS(5)), ' maxtapchange=0'];
Text.Command = ['new regcontrol.Reg6  transformer=Reg6 winding=2  vreg=122  band=2  ptratio=20 ctprim=700  R=3   X=9 TapNum=', num2str(TAPPOS(6)), ' maxtapchange=0'];

Solution.Solve;

if Solution.Converged ~= 1
    disp('Power flow did not converge.')
end

%% Inspecting every nodes to extract phase voltages
ArrangedVoltages_Time = zeros(ActNumBuses,6);
for ii = 1: ActNumBuses

    if sum(phasesinbuses(ii,:)) == 3
        Circuit.SetActiveBus(ArrangedBusName{ii,1});      % It does not matter; choose first element
        ArrangedVoltages_Time(ii,:) = Bus.puVoltages;
    elseif sum(phasesinbuses(ii,:)) == 2
        if phasesinbuses(ii,1) == 1 && phasesinbuses(ii,2) == 1
            Circuit.SetActiveBus(ArrangedBusName{ii,1});  % Choose first/second element
            if numel(Bus.puVoltages) == 4
                ArrangedVoltages_Time(ii,1:4) = Bus.puVoltages;
            else
                ArrangedVoltages_Time(ii,1:4) = Bus.puVoltages(1:4);
            end
        elseif phasesinbuses(ii,1) == 1 && phasesinbuses(ii,3) == 1
            Circuit.SetActiveBus(ArrangedBusName{ii,1});  % Choose first/thirs element
            if numel(Bus.puVoltages) == 4
                ArrangedVoltages_Time(ii,1:2) = Bus.puVoltages(1:2);
                ArrangedVoltages_Time(ii,5:6) = Bus.puVoltages(3:4);
            else
                ArrangedVoltages_Time(ii,1:2) = Bus.puVoltages(1:2);
                ArrangedVoltages_Time(ii,5:6) = Bus.puVoltages(5:6);
            end
        elseif phasesinbuses(ii,2) == 1 && phasesinbuses(ii,3) == 1
            Circuit.SetActiveBus(ArrangedBusName{ii,2});  % Choose second/third element
            if numel(Bus.puVoltages) == 4
                ArrangedVoltages_Time(ii,3:6) = Bus.puVoltages;
            else
                ArrangedVoltages_Time(ii,3:6) = Bus.puVoltages(3:6);
            end
                
        end
    elseif sum(phasesinbuses(ii,:)) == 1
        if phasesinbuses(ii,1) == 1
            Circuit.SetActiveBus(ArrangedBusName{ii,1});  % Choose first element
            ArrangedVoltages_Time(ii,1:2) = Bus.puVoltages(1:2);
        elseif phasesinbuses(ii,2) == 1
            Circuit.SetActiveBus(ArrangedBusName{ii,2});  % Choose first element
            if numel(Bus.puVoltages) == 2
                ArrangedVoltages_Time(ii,3:4) = Bus.puVoltages;
            else
                ArrangedVoltages_Time(ii,3:4) = Bus.puVoltages(3:4);
            end
        else
            Circuit.SetActiveBus(ArrangedBusName{ii,3});  % Choose first element
            if numel(Bus.puVoltages) == 2
                ArrangedVoltages_Time(ii,5:6) = Bus.puVoltages;
            elseif numel(Bus.puVoltages) == 4
                ArrangedVoltages_Time(ii,5:6) = Bus.puVoltages(3:4);
            else
                ArrangedVoltages_Time(ii,5:6) = Bus.puVoltages(5:6);
            end
        end
    end
end

% PU votages
ArrangedVoltagesA = sqrt(sum(ArrangedVoltages_Time(:,1:2).^2,2)); 
ArrangedVoltagesB = sqrt(sum(ArrangedVoltages_Time(:,3:4).^2,2)); 
ArrangedVoltagesC = sqrt(sum(ArrangedVoltages_Time(:,5:6).^2,2)); 

Text.Command = 'Clear';
end
