function [C, Ceq] = nonl_q(x,sys,io,settings)
%Constraint function of torque loop
warning off

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

% The open loop system %
sysol = -Controller*sysdol;

[Gm,Pm,~,Wpm] = margin(sysol);

pm_base = settings.PMLimit;
wpm_base = settings.BandWidth;
gm_base =settings.GMLimit;

if ~isfinite(Pm)|| isnan(Pm)
    Pm = 0;
end

if ~isfinite(Wpm)|| isnan(Wpm)
    Wpm = 0;
end

if ~isfinite(Gm)|| isnan(Gm)
    Gm = 0.1;
end


C(1) = pm_base - Pm;
C(2) = wpm_base - Wpm;
C(3) = gm_base - 20 * log10(Gm);

Ceq = [];


    





