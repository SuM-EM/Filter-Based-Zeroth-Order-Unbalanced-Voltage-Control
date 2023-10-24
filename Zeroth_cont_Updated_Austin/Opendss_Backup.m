clc; clear all;

%% Run at the beginning of the session
DSS = DSS_MATLAB.IDSS;

%% Used aliases
Text = DSS.Text;
Circuit = DSS.ActiveCircuit;
Solution = DSS.ActiveCircuit.Solution;
Bus = DSS.ActiveCircuit.ActiveBus;
Load = DSS.ActiveCircuit.Loads;

%% Clear command interface
Text.Command = 'Clear';
Text.Command = 'Redirect [C:\Users\SuMLiGhT\TAMU\Austin-Grid_Implementation\syn-Austin-TDgrid-v03\syn-austin-D_only-v03\P5U\base\opendss\p5uhs0_1247\p5uhs0_1247--p5udt119\Master.dss]'; 

% Text.Command = ['new Generator.GEN_1 bus=p5ulv138 kw=0 kvar=0']; 
% Obj.ActiveCircuit.Solution.SolveAll
% 
% Obj.Circuits.AllNodeVmagByPhase

if Solution.Converged == 1
    disp('Power flow converged.')
else
    disp('Power flow did not converge.')
end


%% Inspecting every nodes to extract phase voltages
Allbusname = Circuit.AllBusNames;
Allbusbyphases = Circuit.AllNodeNames; numberofphases = numel(Allbusbyphases);
ArrangedBusName = cell(numel(Allbusname),3);

for ii = 1: numberofphases
    str = Allbusbyphases{ii};
    Extract = split(str,'.'); A = Extract{1}; BB = str2double(Extract{2});
    if (ii > 1 && ii <=10032) || ii >= 11039
        if strcmp(A,Aold) == 1
            ArrangedBusName{AA,BB} = str;
        else
            AA = AA + 1; 
            ArrangedBusName{AA,BB} = str;
        end
    elseif ii == 1
        AA = 1;
        ArrangedBusName{AA,BB} = str;
    
    end
    Aold = A;
end

phasesinbuses = zeros(numel(Allbusname),3);
for ii = 1: numel(Allbusname)
    phasesinbuses(ii,:) = 1-cellfun(@isempty,ArrangedBusName(ii,:));
end



ArrangedVoltages = zeros(numel(Allbusname),6);
for ii = 1: numel(Allbusname)

    if sum(phasesinbuses(ii,:)) == 3
        Circuit.SetActiveBus(ArrangedBusName{ii,1});      % It does not matter; choose first element
        ArrangedVoltages(ii,:) = Bus.puVoltages;
    elseif sum(phasesinbuses(ii,:)) == 2
        if phasesinbuses(ii,1) == 1 && phasesinbuses(ii,2) == 1
            Circuit.SetActiveBus(ArrangedBusName{ii,1});  % Choose first/second element
            ArrangedVoltages(ii,1:4) = Bus.puVoltages;
        elseif phasesinbuses(ii,1) == 1 && phasesinbuses(ii,3) == 1
            Circuit.SetActiveBus(ArrangedBusName{ii,1});  % Choose first/thirs element
            ArrangedVoltages(ii,1:2) = Bus.puVoltages(1:2);
            ArrangedVoltages(ii,5:6) = Bus.puVoltages(3:4);
        elseif phasesinbuses(ii,2) == 1 && phasesinbuses(ii,3) == 1
            Circuit.SetActiveBus(ArrangedBusName{ii,2});  % Choose second/third element
            ArrangedVoltages(ii,3:6) = Bus.puVoltages;
        end
    elseif sum(phasesinbuses(ii,:)) == 1
        if phasesinbuses(ii,1) == 1
            Circuit.SetActiveBus(ArrangedBusName{ii,1});  % Choose first element
            ArrangedVoltages(ii,1:2) = Bus.puVoltages;
        elseif phasesinbuses(ii,2) == 1
            Circuit.SetActiveBus(ArrangedBusName{ii,2});  % Choose first element
            ArrangedVoltages(ii,3:4) = Bus.puVoltages;
        else
            Circuit.SetActiveBus(ArrangedBusName{ii,3});  % Choose first element
            ArrangedVoltages(ii,5:6) = Bus.puVoltages;
        end
    end
end

% % Plot A-Phase Voltages
% 
% ArrangedVoltagesPhA = sqrt(sum(ArrangedVoltages(:,1:2).^2,2)); ArrangedVoltagesPhA(ArrangedVoltagesPhA==0) = [];
% 
% figure;
% plot(1:numel(ArrangedVoltagesPhA), ArrangedVoltagesPhA, '-')
% xlim([1 numel(ArrangedVoltagesPhA)])
% 
% ylabel('Voltage (pu)');

% %% Plot B-Phase Voltages
% 
% ArrangedVoltagesPhB = sqrt(sum(ArrangedVoltages(:,3:4).^2,2)); ArrangedVoltagesPhB(ArrangedVoltagesPhB==0) = [];
% 
% figure;
% plot(1:numel(ArrangedVoltagesPhB), ArrangedVoltagesPhB, '-')
% xlim([1 numel(ArrangedVoltagesPhB)])
% 
% ylabel('Voltage (pu)');
% 
% %% Plot C-Phase Voltages
% 
% ArrangedVoltagesPhC = sqrt(sum(ArrangedVoltages(:,5:6).^2,2)); ArrangedVoltagesPhC(ArrangedVoltagesPhC==0) = [];
% 
% figure;
% plot(1:numel(ArrangedVoltagesPhC), ArrangedVoltagesPhC, '-')
% xlim([1 numel(ArrangedVoltagesPhC)])
% 
% ylabel('Voltage (pu)');

%% Select PV locations in Phase A

location = randi([0 1],numel(Allbusname),1); % Randomly select PV locations
capacity = 50*rand(numel(Allbusname),1);     % Randomly select PV capacity
PVinPhaseA = phasesinbuses(:,1).*location;   PVinPhaseA(1) = 0;
PVCapacityinPhaseA = PVinPhaseA.*capacity;

PVData = [0    0    0     0     0    0   0.1    0.2     0.3     0.5     0.8    1.0     0.9    0.8   0.7   0.5   0.5   0.4   0.2   0    0    0    0   0];

for tt = 1: numel(PVData)

DSS.ClearAll();
Text.Command = 'Redirect [C:\Users\SuMLiGhT\TAMU\Austin-Grid_Implementation\syn-Austin-TDgrid-v03\syn-austin-D_only-v03\P5U\base\opendss\p5uhs0_1247\p5uhs0_1247--p5udt119\Master.dss]';

Text.Command=['set mode =snap loadmult=' num2str(1)];   %% Control overall loads here

PV_num = 1;
for ii = 1: numel(Allbusname)
    if PVinPhaseA(ii) == 1
        Text.Command = ['new Generator.PV_A_', num2str(PV_num), ' bus=', ArrangedBusName{ii,1}, ' kw=', num2str(PVCapacityinPhaseA(ii)*PVData(tt)), ' kvar=0'];
        PV_num = PV_num + 1;
    end
end

Solution.Solve;

if Solution.Converged == 1
    disp('Power flow converged.')
else
    disp('Power flow did not converge.')
end


% Store Phase A voltage magnitude

ArrangedVoltages_Time = zeros(numel(Allbusname),6);
for ii = 1: numel(Allbusname)

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
            if numel(Bus.puVoltages) == 2
                ArrangedVoltages_Time(ii,1:2) = Bus.puVoltages;
            else
                ArrangedVoltages_Time(ii,1:2) = Bus.puVoltages(1:2);
            end
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
            else
                ArrangedVoltages_Time(ii,5:6) = Bus.puVoltages(5:6);
            end
        end
    end
end

ArrangedVoltagesA_Time = sqrt(sum(ArrangedVoltages_Time(:,1:2).^2,2)); ArrangedVoltagesA_Time(ArrangedVoltagesA_Time==0) = [];
ArrangedVoltagesPhA_Time(:,tt) = ArrangedVoltagesA_Time;

ArrangedVoltagesB_Time = sqrt(sum(ArrangedVoltages_Time(:,3:4).^2,2)); ArrangedVoltagesB_Time(ArrangedVoltagesB_Time==0) = [];
ArrangedVoltagesPhB_Time(:,tt) = ArrangedVoltagesB_Time;

ArrangedVoltagesC_Time = sqrt(sum(ArrangedVoltages_Time(:,5:6).^2,2)); ArrangedVoltagesC_Time(ArrangedVoltagesC_Time==0) = [];
ArrangedVoltagesPhC_Time(:,tt) = ArrangedVoltagesC_Time;

end
