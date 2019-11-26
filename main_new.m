close all
clear all
clc

warning off

diary on

mkdir('.','PIDCal');

if exist('.\PIDCal\log.txt','file')
    delete('.\PIDCal\log.txt')
end

diary('.\PIDCal\log.txt')

%% ================= File Check ==================== %%
msgID_1 = 'FileCheckError:';

% Required File Check %
msgID_1_1 = [msgID_1, 'FileNotFound'];

if ~exist('.\LinearModel\linmod1.mat','file')
    errmsg_1_1_1 = [datestr(now),'Error:Missing mat file, please check'];
    error(msgID_1_1,errmsg_1_1_1);
end

if ~exist('Settings.xlsx','file')
    errmsg_1_1_2 = [datestr(now), 'Error:Missing setting file, please check'];
    error(msgID_1_1,errmsg_1_1_2);
end

disp([datestr(now),':  初始文件检查正常']);

% Varible Array Check %
msgID_1_2 = [msgID_1 , 'VariableNotFound'];

load('.\LinearModel\linmod1.mat');

v_1_ok = exist('Windspeeds', 'var');
v_2_ok = exist('SYSTURB', 'var');
v_3_ok = exist('PitchAngles', 'var');
v_4_ok = exist('Gbx', 'var');

if ~(v_1_ok&&v_2_ok&&v_3_ok&&v_4_ok)   
    errmsg_1_2 = [datestr(now),'Error:Pivatal variable(s) may missing in the mat file(possibly windspeed, Systurb, pitchangles or gbx), please upload the right file'];
    error(msgID_1_2,errmsg_1_2);
end

disp([datestr(now),':  变量检查正常']);

% Output Check %
msgID_1_3 = [msgID_1, 'OutputNotFound'];

num_output = size(SYSTURB.outputname);
out_pos_Q = 0;
out_pos_P = 0;
out_pos_GS = 0;

for i=1:num_output(1) 
    if strcmpi(strtrim(SYSTURB.outputname(i,:)), 'Generator speed') || strcmpi(strtrim(SYSTURB.outputname(i,:)), 'Measured generator speed')
        out_pos_GS = i;
        break
    end
end

if out_pos_GS ==0
    errmsg_1_3_1 = [datestr(now),'Error:Generator speed is not included in the mat outputs, please upload the right file'];
    error(msgID_1_3,errmsg_1_3_1);
end

for i=1:num_output(1) 
    if strcmpi(strtrim(SYSTURB.outputname(i,:)), 'Blade 1 pitch angle')
        out_pos_P = i;
        break
    end
end

if out_pos_P ==0
    errmsg_1_3_2 = [datestr(now),'Error:Blade 1 pitch angle is not included in the mat outputs, please upload the right file'];
    error(msgID_1_3,errmsg_1_3_2);
end
for i=1:num_output(1) 
    if strcmpi(strtrim(SYSTURB.outputname(i,:)), 'Generator torque')|| strcmpi(strtrim(SYSTURB.outputname(i,:)), 'Demanded generator torque')
        out_pos_Q = i;
        break
    end
end

if out_pos_Q ==0
    errmsg_1_3_3 = [datestr(now),'Error:Generator torque is not included in the mat outputs, please upload the right file'];
    error(msgID_1_3,errmsg_1_3_3);
end

disp([datestr(now),':  模型输出检查正常']);

% Input Check %
msgID_1_4 = [msgID_1, 'InputNotFound'];

num_input = size(SYSTURB.inputname);
in_pos_Q = 0;
in_pos_P = 0;
in_pos_WS = 0;

for i=1:num_input(1) 
    if strcmpi(strtrim(SYSTURB.inputname(i,:)), 'Collective wind speed')
        in_pos_WS = i;
        break
    end
end

if in_pos_WS == 0
    errmsg_1_4_1 = [datestr(now),'Error:Collective wind speed is not included in the mat inputs, please upload the right file'];
    error(msgID_1_4,errmsg_1_4_1);
end

for i=1:num_input(1) 
    if strcmpi(strtrim(SYSTURB.inputname(i,:)), 'Collective pitch angle demand')
        in_pos_P = i;
        break
    end
end

if in_pos_P == 0
    errmsg_1_4_2 = [datestr(now),'Error:Collective pitch angle demand is not included in the mat inputs, please upload the right file'];
    error(msgID_1_4,errmsg_1_4_2);
end
for i=1:num_input(1) 
    if strcmpi(strtrim(SYSTURB.inputname(i,:)), 'Collective generator torque demand')|| strcmpi(strtrim(SYSTURB.inputname(i,:)), 'Generator torque demand')
        in_pos_Q = i;
        break
    end
end

if in_pos_Q == 0
    errmsg_1_4_3 = [datestr(now),'Error:Collective generator torque demand is not included in the mat inputs, please upload the right file'];
    error(msgID_1_4,errmsg_1_4_3);
end

disp([datestr(now),':  模型输入检查正常']);


%% ================= Prepare for tuning ==================== %%
% Init Filters
QFilter = 1;
PFilter = 1;
NFilter = 1;
  
% Time Steps
TsQ = 0.02;	%Torque controller timestep
TsP = 0.02;	%Pitch controller timestep
TsC = 0.02; %Controller platform timestep

% Delays
DelayW = 0.02;	%Speed input delay
DelayA = 0.02;	%Nacelle-x acceleration input delay
DelayL = 0.02;	%Load measurements input delay
DelayQ = 0.06;	%Torque output delay
DelayP = 0.08;	%Pitch output delay

% Define controller properties
CP_type.TimeStep = TsQ;
CP_type.InputDelay = DelayW;
CP_type.OutputDelay = TsC + DelayQ;
CP_type.Name = 'Constant power term';
CP_type.Prefix = 'CP';

QPID_type.TimeStep = TsQ;
QPID_type.InputDelay = DelayW;
QPID_type.OutputDelay = TsC + DelayQ;
QPID_type.Name = 'Torque-speed PID';
QPID_type.Prefix = 'QPID';

PPID_type.TimeStep = TsP;
PPID_type.InputDelay = DelayW;
PPID_type.OutputDelay = TsC + DelayP;
PPID_type.Name = 'Pitch-speed PID';
PPID_type.Prefix = 'PPID';

NAF_type.TimeStep = TsP;
NAF_type.InputDelay = DelayA;
NAF_type.OutputDelay = TsC + DelayP;
NAF_type.Name = 'Nacelle Acceleration Feedback';
NAF_type.Prefix = 'NAF';


%% ================= Load Config ==================== %%
[~, ~, config] = xlsread('Settings.xlsx', 'Config');
msgID_2 = 'SettingError:';
msgID_2_2 = [msgID_2 , 'ValueError'];

% Pitch loop setting %
% Kp setting %
settings.Kpmin = getValue('Kpmin', config);
settings.Kpmax = getValue('Kpmax', config);

if settings.Kpmax < settings.Kpmin
    errmsg_2_2_1 = [datestr(now),'Error:Upper bound should not be smaller than the lower bound for Kp'];
    error(msgID_2_2,errmsg_2_2_1);
elseif settings.Kpmax < 0 || settings.Kpmin < 0
    errmsg_2_2_2 = [datestr(now),'Error:Illegal setting (negative value) for Kp'];
    error(msgID_2_2,errmsg_2_2_2);
end

% Ti setting %
settings.Timin = getValue('Timin', config);
settings.Timax = getValue('Timax', config);

if settings.Timax < settings.Timin
    errmsg_2_2_3 = [datestr(now),'Error:Upper bound should not be smaller than the lower bound for Ti'];
    error(msgID_2_2,errmsg_2_2_3);
elseif settings.Timax < 0 || settings.Timin < 0
    errmsg_2_2_4 = [datestr(now),'Error:Illegal setting (negative value) for Ti'];
    error(msgID_2_2,errmsg_2_2_4);
end

% Kd setting %
settings.isDerEnable = getValue('isDerEnable', config);
if settings.isDerEnable == 1
    settings.Kdmin = getValue('Kdmin', config);
    settings.Kdmax = getValue('Kdmax', config);
    
    if settings.Kdmax < settings.Kdmin
        errmsg_2_2_5 = [datestr(now),'Error:Upper bound should not be smaller than the lower bound for Kd'];
        error(msgID_2_2,errmsg_2_2_5);
    elseif settings.Kdmax < 0 || settings.Kdmin < 0
        errmsg_2_2_6 = [datestr(now),'Error:Illegal setting (negative value) for Kd'];
        error(msgID_2_2,errmsg_2_2_6);
    end 
else
    settings.Kdmin = 0;
    settings.Kdmax = 0;
end

% Td setting %
settings.Td = getValue('Td', config);
if settings.Td <= 0
    errmsg_2_2_7 = [datestr(now),'Error:Illegal setting (negative value) for Td'];
    error(msgID_2_2,errmsg_2_2_7);
end

% NAF loop Setting %
settings.isNAFEnable = getValue('isNAFEnable', config);
settings.NAF_Gain = getValue('NAFGain', config);

if settings.isNAFEnable ~= 1 && settings.isNAFEnable ~= 0
    errmsg_2_2_8 = [datestr(now),'Error:Illegal setting for isNAFEnable'];
    error(msgID_2_2,errmsg_2_2_8);
elseif settings.NAF_Gain <= 0
    errmsg_2_2_9 = [datestr(now),'Error:Illegal setting (negative value) for NAF_Gain'];
    error(msgID_2_2,errmsg_2_2_9);
end

if settings.isNAFEnable == 1
    NAF = tf(1,[1 0]);
    NAF = settings.NAF_Gain*ss(NAF);
    NAF = NAF*NFilter;

    out_pos_NAF = 0;
    for i=1:num_output(1) 
        if strcmpi(strtrim(SYSTURB.outputname(i,:)), 'Nacelle fore-aft acceleration')
            out_pos_NAF = i;
            break
        end
    end

    if out_pos_NAF ==0
        errmsg_1_3_4 = [datestr(now),'Error:Nacelle fore-aft acceleration is not included in the mat outputs, please upload the right file'];
        error(msgID_1_3,errmsg_1_3_4);
    end

    NAF.inputname = strtrim(SYSTURB.outputname(out_pos_NAF,:));
    NAF.outputname = strtrim(SYSTURB.inputname(in_pos_P,:));
end

% Torque loop setting %
% Kp setting %
settings.Qkpmin = getValue('Qkpmin', config);
settings.Qkpmax = getValue('Qkpmax', config);

if settings.Qkpmax < settings.Qkpmin
    errmsg_2_2_10 = [datestr(now),'Error:Upper bound should not be smaller than the lower bound for Qkp'];
    error(msgID_2_2,errmsg_2_2_10);
elseif settings.Qkpmax < 0 || settings.Qkpmin < 0
    errmsg_2_2_11 = [datestr(now),'Error:Illegal setting (negative value) for Qkp'];
    error(msgID_2_2,errmsg_2_2_11);
end

% Ti setting %
settings.Qtimin = getValue('Qtimin', config);
settings.Qtimax = getValue('Qtimax', config);

if settings.Qtimax < settings.Qtimin
    errmsg_2_2_12 = [datestr(now),'Error:Upper bound should not be smaller than the lower bound for Qti'];
    error(msgID_2_2,errmsg_2_2_12);
elseif settings.Qtimax < 0 || settings.Qtimin < 0
    errmsg_2_2_13 = [datestr(now),'Error:Illegal setting (negative value) for Qti'];
    error(msgID_2_2,errmsg_2_2_13);
end

% Margin Limit Setting %
settings.PMLimit = getValue('PMLimit', config);
settings.GMLimit = getValue('GMLimit', config);

if settings.PMLimit <= 0 || settings.GMLimit <=0
    errmsg_2_2_14 = [datestr(now),'Error:Illegal setting (negative value) for margin limits'];
    error(msgID_2_2,errmsg_2_2_14);
end

settings.BandWidth = getValue('BandWidth',config);

if settings.BandWidth <= 0
    errmsg_2_2_15 = [datestr(now),'Error:Illegal setting (negative value) for bandwidth setting'];
    error(msgID_2_2,errmsg_2_2_15);
end

% Optimization Setting %
settings.Gen = getValue('Gen', config);
settings.Pop = getValue('Pop', config);

if mod(settings.Gen,4) ~= 0 || mod(settings.Pop,4) ~= 0
    errmsg_2_2_16 = [datestr(now),'Error:Upper bound should not be smaller than the lower bound for Qti'];
    error(msgID_2_2,errmsg_2_2_16);
elseif settings.Gen <= 0 || settings.Pop <= 0
    errmsg_2_2_17 = [datestr(now),'Error:Illegal setting (negative value) for optimization setting'];
    error(msgID_2_2,errmsg_2_2_17);
end

settings.FilterSet = getValue('FilterSet', config);

if settings.FilterSet == 1
    filter_file = 'Settings.xlsx';
elseif settings.FilterSet == 0
    filter_file = '\Temp\Filters.xlsx';
else
    errmsg_2_2_18 = [datestr(now),'Error:Illegal setting (negative value) for filter setting'];
    error(msgID_2_2,errmsg_2_2_18);
end

settings.doComfit = getValue('doComfit', config);
settings.comfitindex = getValue('comfitindex', config);

if settings.doComfit ~= 1 && settings.doComfit ~= 0
    errmsg_2_2_19 = [datestr(now),'Error:Illegal setting (negative value) for comfit setting'];
    error(msgID_2_2,errmsg_2_2_19);
elseif settings.comfitindex <= 0 
    errmsg_2_2_20 = [datestr(now),'Error:Illegal setting (negative value) for comfit setting'];
    error(msgID_2_2,errmsg_2_2_20);
end

settings.TorqueOnly = getValue('TorqueOnly', config);
settings.PitchOnly = getValue('PitchOnly', config);

if settings.TorqueOnly ~= 1 && settings.TorqueOnly ~= 0
    errmsg_2_2_21 = [datestr(now),'Error:Illegal setting (negative value) for torque loop only setting'];
    error(msgID_2_2,errmsg_2_2_21);
elseif settings.PitchOnly ~= 1 && settings.PitchOnly ~= 0
    errmsg_2_2_22 = [datestr(now),'Error:Illegal setting (negative value) for pitch loop only setting'];
    error(msgID_2_2,errmsg_2_2_22);
end

disp([datestr(now),':  配置参数检查正常']);

%% ================= Station Setting ==================== %%

Order = 1;  %Station Number
S_sign = 1;
FinePitch = PitchAngles(1); %Minimum Fine PitchAngles
PitchLoop = 0; %PitchAngles Loop Control Flag

windgap = Windspeeds(2) - Windspeeds(1);

if windgap <= 0.5
    StationStep = 4;%Gap between two stations
elseif windgap >= 2
    StationStep = 1;
else
    StationStep = 2;
end

while ~PitchLoop&&S_sign < length(PitchAngles)
    if PitchAngles(S_sign) > FinePitch
        Q_station = S_sign - 1; %% Torque loop
        Q_pitchangle = PitchAngles(S_sign - 1);
        Q_windspeed = Windspeeds(S_sign - 1);
        
        P_station(Order) = S_sign; %%Pitch loop first station
        P_pitchangle(Order) = PitchAngles(S_sign);
        P_windspeed(Order) = Windspeeds(S_sign);
        
        PitchLoop = 1;
    end
    S_sign = S_sign + 1;
end

if StationStep == 4
    Order = Order + 1;
    P_station(Order) = S_sign + 1;
    P_pitchangle(Order) = PitchAngles(S_sign + 1);  %%Second Station
    P_windspeed(Order) = Windspeeds(S_sign + 1);
end
    
for j = S_sign + StationStep - 1:StationStep:length(PitchAngles)
    Order = Order + 1;
    if j <= (length(PitchAngles)-StationStep)
        P_station(Order) = j;
        P_pitchangle(Order) = PitchAngles(j);  %%Further Station
        P_windspeed(Order) = Windspeeds(j);
    else
        P_station(Order) = length(PitchAngles);
        P_pitchangle(Order) = PitchAngles(end);  %%Final Station
        P_windspeed(Order) = Windspeeds(end);
    end
end

mkdir('.\PIDCal','Torque');
station_size = length(P_station);
for i = 1 : station_size   
    mkdir('.\PIDCal', ['Pitch', num2str(i)]);
end

%% ================= Filters Loading ==================== %%
filters_p = xlsread(filter_file, 'Pitch');

w = filters_p(1,3);  z= filters_p(1,4);
num = w^2;
den = [1 2*z*w w^2];
Plpf = tf(num,den);
PFilter = PFilter*ss(Plpf);

for i = 2:length(filters_p) 

    wn = filters_p(i,1); dn = filters_p(i,2); wd = filters_p(i,3); dd = filters_p(i,4);
    num = [1/wn.^2 2*dn/wn 1];
    den = [1/wd.^2 2*dd/wd 1];
    n = ss(tf(num,den));
    PFilter = PFilter*n;

end

filters_q = xlsread(filter_file, 'Torque');

w = filters_q(1,3);  z= filters_q(1,4);
num = w^2;
den = [1 2*z*w w^2];
Qlpf = tf(num,den);
QFilter = QFilter*ss(Qlpf);

for i = 2:length(filters_q) 

    wn = filters_q(i,1); dn = filters_q(i,2); wd = filters_q(i,3); dd = filters_q(i,4);
    num = [1/wn.^2 2*dn/wn 1];
    den = [1/wd.^2 2*dd/wd 1];
    n = ss(tf(num,den));
    QFilter = QFilter*n;

end

filters_naf = xlsread(filter_file, 'NAF');

w = filters_naf(1,3);  z= filters_naf(1,4);
num = w^2;
den = [1 2*z*w w^2];
Nlpf = tf(num,den);
NFilter = NFilter*ss(Nlpf);

for i = 2:length(filters_naf) 

    wn = filters_naf(i,1); dn = filters_naf(i,2); wd = filters_naf(i,3); dd = filters_naf(i,4);
    num = [1/wn.^2 2*dn/wn 1];
    den = [1/wd.^2 2*dd/wd 1];
    n = ss(tf(num,den));
    NFilter = NFilter*n;

end

disp([datestr(now),':  滤波器已加载']);


%% ================= Handle Initializaition ==================== %%
handle = zeros(3, length(P_station)+1);
handle(2,2) = 1;
for i = 1:length(P_station)
    handle(1,i+1) = 1;
end

%% ================= Loop Tuning (Step 1)==================== %%
% Define input and output %
io.q_input = strtrim(SYSTURB.outputname(out_pos_GS,:));
io.q_output = strtrim(SYSTURB.inputname(in_pos_Q,:));
io.q_input_disp = strtrim(SYSTURB.inputname(in_pos_WS,:)); 
io.q_output_disp = strtrim(SYSTURB.outputname(out_pos_Q,:)); 

io.p_input = strtrim(SYSTURB.outputname(out_pos_GS,:));
io.p_output = strtrim(SYSTURB.inputname(in_pos_P,:));
io.p_input_disp = strtrim(SYSTURB.inputname(in_pos_WS,:)); 
io.p_output_disp = strtrim(SYSTURB.outputname(out_pos_P,:)); 

io.PPID_type = PPID_type;
io.QPID_type = QPID_type;
io.PFilter = PFilter;
io.QFilter = QFilter;

% Experiment start %

for i = 1:2
    if handle(1,i) == 0 && ~settings.PitchOnly
        disp([datestr(now),':  扭矩回路整定开始']);
        
        sys = ss(SYSTURB.A(:,:,Q_station,1),SYSTURB.B(:,:,Q_station,1),SYSTURB.C(:,:,Q_station,1),SYSTURB.D(:,:,Q_station,1));
        sys.inputname = SYSTURB.inputname;  
        sys.outputname = SYSTURB.outputname; 
        
        move = 0;
        
        result = GAMultiObj(sys,handle(1,i),handle(2,i),handle(3,i),settings,io);
        
        pack;
        
        
        if ~isempty(find(result(:,3)>90, 1))
            move = 1;
            
            sys = ss(SYSTURB.A(:,:,Q_station+1,1),SYSTURB.B(:,:,Q_station+1,1),SYSTURB.C(:,:,Q_station+1,1),SYSTURB.D(:,:,Q_station+1,1));
            sys.inputname = SYSTURB.inputname;  
            sys.outputname = SYSTURB.outputname; 

            result = GAMultiObj(sys,handle(1,i),handle(2,i),handle(3,i),settings,io);
            
            pack;
        end
        
        r_fre = zeros(size(result,1),4);
        r_score = zeros(size(result,1),1);
        r_id = zeros(size(result,1),1);
        
        if move
            r_station = PitchAngles(Q_station+1)*ones(size(result,1),1);
        else
            r_station = PitchAngles(Q_station)*ones(size(result,1),1);
        end
        
        for j = 1:size(result,1)
            Kp = result(j,1);
            Ti = result(j,2);
            Kd = 0;
            Td = 0;
            
            QPID = ss(tf([Kd + Kp*Td, Kp + Kp/Ti*Td, Kp/Ti], [Td, 1, 0]));
            QPID = QPID*QFilter; 
            QPID.inputname = io.q_input;
            QPID.outputname = io.q_output;

            % Assign Loop %
            Controller = QPID;
            ControllerType = QPID_type;

            sysd = convertSsTime(sys,ControllerType.TimeStep,'ZOH','Plant');
            sysdol = getSiso(sysd,Controller.outputname,Controller.inputname);
            sysdol.outputdelay = round(ControllerType.InputDelay/ControllerType.TimeStep);
            sysdol.inputdelay = round((ControllerType.InputDelay + ControllerType.OutputDelay)/ControllerType.TimeStep)-sysdol.outputdelay;
            Controller = c2d(Controller,ControllerType.TimeStep,'tustin');   

            % The open loop system %
            sysol = -Controller*sysdol;

            [Gm,Pm,Wgm,Wpm] = margin(sysol);
            
            r_fre(j,1) = 20 * log10(Gm);
            r_fre(j,2) = Pm;
            r_fre(j,3) = Wgm;
            r_fre(j,4) = Wpm;
            
            Controller.inputname = QPID.inputname;
            Controller.outputname = QPID.outputname;

            ControllerDelay = Controller;
            ControllerDelay.inputdelay = round(ControllerType.InputDelay/ControllerType.TimeStep); % add delays to the controller
            ControllerDelay.outputdelay = round((ControllerType.InputDelay + ControllerType.OutputDelay)/ControllerType.TimeStep)-sysdol.outputdelay; % add delays to the controller

            % The closed loop system %
            syscl = addController(sys, ControllerDelay, ControllerType);
            
            r_score(j) = grading(result(j,3),result(j,4),result(j,5),result(j,6));
            r_id(j) = j;
            
            h = figure('Visible', 'off');
            step(getSiso(syscl,io.q_input_disp,io.q_output_disp),100);
            xlabel Time
            ylabel GeneratorTorque
            grid on
            saveas(h, ['.\PIDCal\Torque\Torque_', num2str(j),'_step'],'png')
            
            g = figure('Visible', 'off');
            margin(sysol)
            grid on
            saveas(g, ['.\PIDCal\Torque\Torque_', num2str(j),'_bode'],'png')
        end
        
        q_result = [r_id,r_station,result,r_fre,r_score];
        
        filename = '.\PIDCal\Torque\Result.xlsx';

        if exist(filename,'file')
            delete(filename)
        end
        
        Label_q = {'ID','Station','Kp','Ti',...
                   'Settling Time','Rise Time','OverShoot','Oscillation',...
                   'Gain Margin','Phase Margin','Gain Margin Frequency','Phase margin frequency',...
                   'Score'};
               
        xlswrite(filename,Label_q,1);
        xlswrite(filename,q_result,1,'A2');
        
        [c_h,~] = find(r_score == max(r_score));
        
        for j = 1:length(c_h)
            r_b = [100,100];
            r_order = c_h(j);
            
            if (result(r_order,3) + 5*result(r_order,4)) <= (r_b(1) + 5*r_b(2))
                r_b(1) = result(r_order,3);
                r_b(2) = result(r_order,4);
                r_ob = r_order;
            end
        end
        
        q_rb = [r_id(r_ob),r_station(r_ob),result(r_ob,:),r_fre(r_ob,:),r_score(r_ob)];
        
        filename = '.\PIDCal\Torque\Result.txt';
        
        if exist(filename,'file')
            delete(filename)
        end
        
        fid = fopen(filename,'wt');
        q_rf = result(r_ob,:);
        
        for j = 1:size(q_rf,1)
            for k = 1:size(q_rf,2)
                fprintf(fid, '%.4f ', q_rf(j,k));
            end
            if j < size(q_rf,1)
                fprintf(fid, '\n');
            end
        end
        fclose(fid);
        disp([datestr(now),':  扭矩整定完成']);
       
    elseif handle(1,i) == 1 && ~settings.TorqueOnly
        
        clearvars -except SYSTURB settings io handle PitchAngles P_pitchangle P_station i NAF NAF_type Label_q q_rb
        
        disp([datestr(now),':  变桨回路整定开始']);
        
        sys = ss(SYSTURB.A(:,:,P_station(1),1),SYSTURB.B(:,:,P_station(1),1),SYSTURB.C(:,:,P_station(1),1),SYSTURB.D(:,:,P_station(1),1));
        sys.inputname = SYSTURB.inputname;  
        sys.outputname = SYSTURB.outputname;
        
        if settings.isNAFEnable == 1
            sys = addController(sys, NAF, NAF_type);
        end 
        
        result = GAMultiObj(sys,handle(1,i),handle(2,i),handle(3,i),settings,io);
        
        pack;
        
        r_fre = zeros(size(result,1),4);
        r_score = zeros(size(result,1),1);
        r_id = zeros(size(result,1),1);
        r_station = PitchAngles(P_station(1))*ones(size(result,1),1);
        
        for j = 1:size(result,1)
            Kp = result(j,1);
            Ti = result(j,2);
            Kd = result(j,3);
            Td = result(j,4);
            
            PPID = ss(tf([Kd + Kp*Td, Kp + Kp/Ti*Td, Kp/Ti], [Td, 1, 0]));
            PPID = PPID*io.PFilter; 
            PPID.inputname = io.p_input;
            PPID.outputname = io.p_output;

            % Assign Loop %
            Controller = PPID;
            ControllerType = io.PPID_type;

            sysd = convertSsTime(sys,ControllerType.TimeStep,'ZOH','Plant');
            sysdol = getSiso(sysd,Controller.outputname,Controller.inputname);
            sysdol.outputdelay = round(ControllerType.InputDelay/ControllerType.TimeStep);
            sysdol.inputdelay = round((ControllerType.InputDelay + ControllerType.OutputDelay)/ControllerType.TimeStep)-sysdol.outputdelay;
            Controller = c2d(Controller,ControllerType.TimeStep,'tustin');   

            % The open loop system %
            sysol = -Controller*sysdol;

            [Gm,Pm,Wgm,Wpm] = margin(sysol);
            
            r_fre(j,1) = 20 * log10(Gm);
            r_fre(j,2) = Pm;
            r_fre(j,3) = Wgm;
            r_fre(j,4) = Wpm;
            
            Controller.inputname = PPID.inputname;
            Controller.outputname = PPID.outputname;

            ControllerDelay = Controller;
            ControllerDelay.inputdelay = round(ControllerType.InputDelay/ControllerType.TimeStep); % add delays to the controller
            ControllerDelay.outputdelay = round((ControllerType.InputDelay + ControllerType.OutputDelay)/ControllerType.TimeStep)-sysdol.outputdelay; % add delays to the controller

            % The closed loop system %
            syscl = addController(sys, ControllerDelay, ControllerType);
            
            r_score(j) = grading(result(j,5),result(j,6),result(j,7),result(j,8));
            r_id(j) = j;
            
            h = figure('Visible', 'off');
            step(getSiso(syscl,io.p_input_disp,io.p_output_disp),100);
            xlabel Time
            ylabel PitchAngle
            grid on
            saveas(h, ['.\PIDCal\Pitch1\Pitch_1_', num2str(j),'_step'],'png')
            
            g = figure('Visible', 'off');
            margin(sysol)
            grid on
            saveas(g, ['.\PIDCal\Pitch1\Pitch_1_', num2str(j),'_bode'],'png')
        end
        
        p_1_result = [r_id,r_station,result,r_fre,r_score];
        
        filename = '.\PIDCal\Pitch1\Result.xlsx';

        if exist(filename,'file')
            delete(filename)
        end
        
        Label_p = {'ID','Station','Kp','Ti','Kd','Td'...
                   'Settling Time','Rise Time','OverShoot','Oscillation',...
                   'Gain Margin','Phase Margin','Gain Margin Frequency','Phase margin frequency',...
                   'Score'};
               
        xlswrite(filename,Label_p,1);
        xlswrite(filename,p_1_result,1,'A2');
        
        [c_h,~] = find(r_score == max(r_score));
        
        for j = 1:length(c_h)
            r_b = [100,100];
            r_order = c_h(j);
            
            if (result(r_order,5) + 5*result(r_order,6)) <= (r_b(1) + 5*r_b(2))
                r_b(1) = result(r_order,5);
                r_b(2) = result(r_order,6);
                r_ob = r_order;
            end
        end
        p_1_rb = [r_id(r_ob),r_station(r_ob),result(r_ob,:),r_fre(r_ob,:),r_score(r_ob)];
        
        filename = '.\PIDCal\Pitch1\Result.txt';
        
        if exist(filename,'file')
            delete(filename)
        end
        
        fid = fopen(filename,'wt');
        p_1_rf = result(r_order,:);
        
        for j = 1:size(p_1_rf,1)
            for k = 1:size(p_1_rf,2)
                fprintf(fid, '%.4f ', p_1_rf(j,k));
            end
            if j < size(p_1_rf,1)
                fprintf(fid, '\n');
            end
        end
        fclose(fid);
        
        disp([datestr(now),':  变桨1整定完成']);
        
    end
end

%% ================= Loop Tuning (Step 2)==================== %%
if ~settings.TorqueOnly
    p_rb = zeros(size(handle,2)-2,size(p_1_rb,2));

    if settings.doComfit == 1 
        c_s = 0;

        for i = 1:length(PitchAngles)
            if PitchAngles(i)>PitchAngles(1)
                c_s = i;
                break;
            end
        end

        pitch_angles = PitchAngles(c_s:end);
        pitch_angles_dot = zeros(1,length(pitch_angles));

        for i = 1:length(pitch_angles)
            if i == 1
                pitch_angles_dot(i) = pitch_angles(i)-PitchAngles(1);
            else
                pitch_angles_dot(i) = pitch_angles(i)-pitch_angles(i-1);
            end
        end

        P_pitchangle_dot = zeros(1,length(P_pitchangle));

        for i  = 1:length(P_pitchangle)
            ind = find(pitch_angles == P_pitchangle(i));
            P_pitchangle_dot(i) = pitch_angles_dot(ind);
        end

        kp_list = zeros(1,length(P_pitchangle));
        ti_list = zeros(1,length(P_pitchangle));
        kd_list = zeros(1,length(P_pitchangle));

        kp_list(1) = p_1_rb(3);
        ti_list(1) = p_1_rb(4);
        kd_list(1) = p_1_rb(5);

        fit_coe = settings.comfitindex;

        for i  = 2:length(P_pitchangle)
            coe = power(P_pitchangle_dot(i)/P_pitchangle_dot(1),fit_coe);
            kp_list(i) = coe*kp_list(1);
            ti_list(i) = coe*ti_list(1);
            kd_list(i) = coe*kd_list(1);

            Kp = kp_list(i);
            Ti = ti_list(i);
            Kd = kd_list(i);
            Td = p_1_rb(6);

            PPID = ss(tf([Kd + Kp*Td, Kp + Kp/Ti*Td, Kp/Ti], [Td, 1, 0]));
            PPID = PPID*io.PFilter; 
            PPID.inputname = io.p_input;
            PPID.outputname = io.p_output;

            % Assign Loop %
            Station = P_station(i);
            iAzimuth = 1;
            sys = ss(SYSTURB.A(:,:,Station,iAzimuth),SYSTURB.B(:,:,Station,iAzimuth),SYSTURB.C(:,:,Station,iAzimuth),SYSTURB.D(:,:,Station,iAzimuth));
            sys.inputname = SYSTURB.inputname;  %A list of all the inputnames
            sys.outputname = SYSTURB.outputname; %A 

            if settings.isNAFEnable == 1
                sys = addController(sys, NAF, NAF_type);
            end

            Controller = PPID;
            ControllerType = io.PPID_type;

            sysd = convertSsTime(sys,ControllerType.TimeStep,'ZOH','Plant');
            sysdol = getSiso(sysd,Controller.outputname,Controller.inputname);
            sysdol.outputdelay = round(ControllerType.InputDelay/ControllerType.TimeStep);
            sysdol.inputdelay = round((ControllerType.InputDelay + ControllerType.OutputDelay)/ControllerType.TimeStep)-sysdol.outputdelay;
            Controller = c2d(Controller,ControllerType.TimeStep,'tustin');  

            sysol = -Controller*sysdol;
            [Gm,Pm,Wgm,Wpm] = margin(sysol);

            Controller.inputname = PPID.inputname;
            Controller.outputname = PPID.outputname;

            ControllerDelay = Controller;
            ControllerDelay.inputdelay = round(ControllerType.InputDelay/ControllerType.TimeStep); % add delays to the controller
            ControllerDelay.outputdelay = round((ControllerType.InputDelay + ControllerType.OutputDelay)/ControllerType.TimeStep)-sysdol.outputdelay; % add delays to the controller

            % The closed loop system %
            syscl = addController(sys, ControllerDelay, ControllerType);
            [y_sam,t_sam] = step(getSiso(syscl,io.p_input_disp,io.p_output_disp),100);
            Step1 = stepinfo(y_sam,t_sam);

            t_s = Step1.SettlingTime;
            osc_num = 0;
            pos_f = find(t_sam>=t_s,1);
            y_f = y_sam(pos_f);

            for s = 2:length(y_sam)
                if (y_sam(s-1) - y_f)*(y_sam(s) - y_f) < 0
                    osc_num = osc_num + 1;
                end
            end

            score = grading(Step1.SettlingTime,Step1.RiseTime,Step1.Overshoot,osc_num);
            id = 1;

            p_result = [id,PitchAngles(P_station(i)),Kp,Ti,Kd,Td,Step1.SettlingTime,Step1.RiseTime,Step1.Overshoot,osc_num,Gm,Pm,Wgm,Wpm,score];

            filename = ['.\PIDCal\Pitch',num2str(i),'\Result.xlsx'];

            if exist(filename,'file')
                delete(filename)
            end

            xlswrite(filename,Label_p,1);
            xlswrite(filename,p_result,1,'A2');

            p_rb(i-1,:) = p_result;

            h = figure('Visible', 'off');
            step(getSiso(syscl,io.p_input_disp,io.p_output_disp),100);
            xlabel Time
            ylabel PitchAngle
            grid on
            saveas(h, ['.\PIDCal\Pitch', num2str(i),'\Pitch_', num2str(i),'_comfit_step'],'png')

            g = figure('Visible', 'off');
            margin(sysol);
            grid on
            saveas(g, ['.\PIDCal\Pitch', num2str(i),'\Pitch_', num2str(i),'_comfit_bode'],'png')

            result = [Kp,Ti,Kd,Td,Step1.SettlingTime,Step1.RiseTime,Step1.Overshoot,osc_num];

            fid = fopen(['.\PIDCal\Pitch', num2str(i),'\Result.txt'],'wt');
            for k = 1:size(result,1)
                for j = 1:size(result,2)
                    fprintf(fid, '%.4f ', result(k,j));
                end
                if k < size(result,1)
                    fprintf(fid, '\n');
                end
            end
            fclose(fid); 
            disp([datestr(now),':  变桨',num2str(i),'整定完成']);
        end
    else
        for i = 2:size(handle,2)-1
            
            clearvars -except SYSTURB settings io handle PitchAngles P_station i p_1_rb NAF NAF_type Label_p p_rb Label_q q_rb
            
            handle(3,i+1) = p_1_rb(6);

            sys = ss(SYSTURB.A(:,:,P_station(i),1),SYSTURB.B(:,:,P_station(i),1),SYSTURB.C(:,:,P_station(i),1),SYSTURB.D(:,:,P_station(i),1));
            sys.inputname = SYSTURB.inputname;  
            sys.outputname = SYSTURB.outputname;

            if settings.isNAFEnable == 1
                sys = addController(sys, NAF, NAF_type);
            end 

            result = GAMultiObj(sys,handle(1,i+1),handle(2,i+1),handle(3,i+1),settings,io);
            
            pack;

            r_fre = zeros(size(result,1),4);
            r_score = zeros(size(result,1),1);
            r_id = zeros(size(result,1),1);
            r_station = PitchAngles(P_station(i))*ones(size(result,1),1);

            for j = 1:size(result,1)
                Kp = result(j,1);
                Ti = result(j,2);
                Kd = result(j,3);
                Td = result(j,4);

                PPID = ss(tf([Kd + Kp*Td, Kp + Kp/Ti*Td, Kp/Ti], [Td, 1, 0]));
                PPID = PPID*io.PFilter; 
                PPID.inputname = io.p_input;
                PPID.outputname = io.p_output;

                % Assign Loop %
                Controller = PPID;
                ControllerType = io.PPID_type;

                sysd = convertSsTime(sys,ControllerType.TimeStep,'ZOH','Plant');
                sysdol = getSiso(sysd,Controller.outputname,Controller.inputname);
                sysdol.outputdelay = round(ControllerType.InputDelay/ControllerType.TimeStep);
                sysdol.inputdelay = round((ControllerType.InputDelay + ControllerType.OutputDelay)/ControllerType.TimeStep)-sysdol.outputdelay;
                Controller = c2d(Controller,ControllerType.TimeStep,'tustin');   

                % The open loop system %
                sysol = -Controller*sysdol;

                [Gm,Pm,Wgm,Wpm] = margin(sysol);

                r_fre(j,1) = Gm;
                r_fre(j,2) = Pm;
                r_fre(j,3) = Wgm;
                r_fre(j,4) = Wpm;

                Controller.inputname = PPID.inputname;
                Controller.outputname = PPID.outputname;

                ControllerDelay = Controller;
                ControllerDelay.inputdelay = round(ControllerType.InputDelay/ControllerType.TimeStep); % add delays to the controller
                ControllerDelay.outputdelay = round((ControllerType.InputDelay + ControllerType.OutputDelay)/ControllerType.TimeStep)-sysdol.outputdelay; % add delays to the controller

                % The closed loop system %
                syscl = addController(sys, ControllerDelay, ControllerType);

                r_score(j) = grading(result(j,5),result(j,6),result(j,7),result(j,8));
                r_id(j) = j;

                h = figure('Visible', 'off');
                step(getSiso(syscl,io.p_input_disp,io.p_output_disp),100);
                xlabel Time
                ylabel PitchAngle
                grid on
                saveas(h, ['.\PIDCal\Pitch', num2str(i),'\Pitch_', num2str(i),'_',num2str(j),'_step'],'png')

                g = figure('Visible', 'off');
                margin(sysol);
                grid on
                saveas(g, ['.\PIDCal\Pitch', num2str(i),'\Pitch_', num2str(i),'_',num2str(j),'_bode'],'png')
            end

            p_result = [r_id,r_station,result,r_fre,r_score];

            filename = ['.\PIDCal\Pitch', num2str(i),'\Result.xlsx'];

            if exist(filename,'file')
                delete(filename)
            end

            xlswrite(filename,Label_p,1);
            xlswrite(filename,p_result,1,'A2');

            [c_h,~] = find(r_score == max(r_score));

            for j = 1:length(c_h)
                r_b = [100,100];
                r_order = c_h(j);

                if (result(r_order,5) + 5*result(r_order,6)) <= (r_b(1) + 5*r_b(2))
                    r_b(1) = result(r_order,5);
                    r_b(2) = result(r_order,6);
                    r_ob = r_order;
                end
            end
            p_rb(i-1,:) = [r_id(r_ob),r_station(r_ob),result(r_ob,:),r_fre(r_ob,:),r_score(r_ob)];

            fid = fopen(['.\PIDCal\Pitch', num2str(i),'\Result.txt'],'wt');
            p_rf = result(r_ob,:);

            for k = 1:size(p_rf,1)
                for j = 1:size(p_rf,2)
                    fprintf(fid, '%.4f ', p_rf(k,j));
                end
                if k < size(p_rf,1)
                    fprintf(fid, '\n');
                end
            end
            fclose(fid);
            disp([datestr(now),':  变桨',num2str(i),'整定完成']);
        end

    end
end

%% ================= Result Collection ==================== %%

filename = '.\PIDCal\Result_all.xlsx';

if ~settings.PitchOnly && ~settings.TorqueOnly
    if exist(filename,'file')
        delete(filename)
    end
end
if ~settings.TorqueOnly
    xlswrite(filename,Label_p,1);
    xlswrite(filename,p_1_rb,1,'A2');
    xlswrite(filename,p_rb,1,'A3');
end
if ~settings.PitchOnly
    xlswrite(filename,Label_q,2);
    xlswrite(filename,q_rb,2,'A2');
end

diary off


