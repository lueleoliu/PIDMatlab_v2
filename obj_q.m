function y = obj_q(x,sys,io)
%Objective function of torque loop
Kp = x(1);
Ti = x(2);
Kd = 0;
Td = 0;

QPID = ss(tf([Kd + Kp*Td, Kp + Kp/Ti*Td, Kp/Ti], [Td, 1, 0]));
QPID = QPID*io.QFilter; 
QPID.inputname = io.q_input;
QPID.outputname = io.q_output;

% Assign Loop %

Controller = QPID;
ControllerType = io.QPID_type;

sysd = convertSsTime(sys,ControllerType.TimeStep,'ZOH','Plant');
sysdol = getSiso(sysd,Controller.outputname,Controller.inputname);
sysdol.outputdelay = round(ControllerType.InputDelay/ControllerType.TimeStep);
sysdol.inputdelay = round((ControllerType.InputDelay + ControllerType.OutputDelay)/ControllerType.TimeStep)-sysdol.outputdelay;
Controller = c2d(Controller,ControllerType.TimeStep,'tustin');   

Controller.inputname = QPID.inputname;
Controller.outputname = QPID.outputname;

ControllerDelay = Controller;
ControllerDelay.inputdelay = round(ControllerType.InputDelay/ControllerType.TimeStep); % add delays to the controller
ControllerDelay.outputdelay = round((ControllerType.InputDelay + ControllerType.OutputDelay)/ControllerType.TimeStep)-sysdol.outputdelay; % add delays to the controller

% The closed loop system %
syscl = addController(sys, ControllerDelay, ControllerType);
[y_sam,t_sam] = step(getSiso(syscl,io.q_input_disp,io.q_output_disp),100);
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
    





