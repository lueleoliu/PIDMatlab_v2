function y = obj_p(x,sys,io)
%Objective function of pitch loop
Kp = x(1);
Ti = x(2);
Kd = x(3);
Td = x(4);

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
y(1) = t_s;
y(2) = Step1.RiseTime;
y(3) = Step1.Overshoot;

osc_num = 0;
pos_f = find(t_sam>=t_s,1);
y_f = y_sam(pos_f);

for s = 2:length(y_sam)
    if (y_sam(s-1) - y_f)*(y_sam(s) - y_f) < 0
        osc_num = osc_num + 1;
    end
end
y(4) = osc_num;
    





