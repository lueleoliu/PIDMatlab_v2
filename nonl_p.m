function [C, Ceq] = nonl_p(x,sys,io,settings)
%Constraint function of pitch loop
warning off

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


    





