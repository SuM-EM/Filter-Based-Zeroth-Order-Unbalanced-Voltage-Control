function [ArrangedBusName,phasesinbuses,ActNumBuses] = opendss_startup(masterfile,soucepath)
%opendss_startup Identify active nodes in the distribution system 
%                Input:  masterfile:      Missing string for master.dss
%                        soucepath:       Path of source folder
%                Output: ArrangedBusName: which phases are present in each buses; inactive buses are removed (valid node numbers)
%                        phasesinbuses:   which phases are present in each buses; inactive buses are removed (identifier)
%                        ActNumBuses:     How many number of buses are active

%% Run at the beginning of the session
DSS = DSS_MATLAB.IDSS;

%% Used aliases
Text = DSS.Text;
Circuit = DSS.ActiveCircuit;
Solution = DSS.ActiveCircuit.Solution;
Bus = DSS.ActiveCircuit.ActiveBus;
Load = DSS.ActiveCircuit.Loads;

%% Clear command interface
DSS.ClearAll();
Text.Command = ['Redirect [',soucepath,masterfile,']']; 

if Solution.Converged ~= 1
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

ArrangedBusName = ArrangedBusName(1:AA,:);   %% Edit arranged bus names; not all bus names are used in the final list
ActNumBuses = AA;                            %% Actual number of buses

phasesinbuses = zeros(ActNumBuses,3);
for ii = 1: ActNumBuses
    phasesinbuses(ii,:) = 1-cellfun(@isempty,ArrangedBusName(ii,:));
end

%% All buses active nodes are loaded with 1kW of PV
%% in Phase A
PVinPhaseA = phasesinbuses(:,1);   PVinPhaseA(1) = 0;
%% in Phase B
PVinPhaseB = phasesinbuses(:,2);   PVinPhaseB(1) = 0;
%% in Phase C
PVinPhaseC = phasesinbuses(:,3);   PVinPhaseC(1) = 0;

Text.Command = 'Clear';
Text.Command = ['Redirect [',soucepath,masterfile,']']; 
%% Define all the PV outputs
PV_num = 1;
for ii = 1: ActNumBuses
    if PVinPhaseA(ii) > 0
        Text.Command = ['new Generator.PV_A_', num2str(PV_num), ' bus=', ArrangedBusName{ii,1}, ' Model=1 kw=1 kvar=1'];
        PV_num = PV_num + 1;
    end
end

PV_num = 1;
for ii = 1: ActNumBuses
    if PVinPhaseB(ii) > 0
        Text.Command = ['new Generator.PV_B_', num2str(PV_num), ' bus=', ArrangedBusName{ii,2}, ' Model=1 kw=1 kvar=1'];
        PV_num = PV_num + 1;
    end
end

PV_num = 1;
for ii = 1: ActNumBuses
    if PVinPhaseC(ii) > 0
        Text.Command = ['new Generator.PV_C_', num2str(PV_num), ' bus=', ArrangedBusName{ii,3}, ' Model=1 kw=1 kvar=1'];
        PV_num = PV_num + 1;
    end
end

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
ArrangedVoltagespu = [ sqrt(sum(ArrangedVoltages_Time(:,1:2).^2,2)), ... 
                            sqrt(sum(ArrangedVoltages_Time(:,3:4).^2,2)), ... 
                                sqrt(sum(ArrangedVoltages_Time(:,5:6).^2,2))]; 
                            

phasesinbuses(ArrangedVoltagespu < 0.01) = 0;
phasesinbuses(ArrangedVoltagespu > 10) = 0;

for ii = 1:ActNumBuses
    for jj = 1:3
        if phasesinbuses(ii,jj) == 0
            ArrangedBusName{ii,jj} = {};
        end
    end
end

Text.Command = 'Clear';
end

