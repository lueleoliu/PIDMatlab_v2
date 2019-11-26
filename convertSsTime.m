function sysOut = convertSsTime(sysIn,timeStep,method,sysName)
%Converts a state space variable to a different time domain.
%   
%   SYSOUT = CONVERTSSTIME(SYSIN,TIMESTEP,METHOD,SYSNAME) uses C2D, D2C  
%   or D2D as appropriate (with the 'METHOD' option) to convert the state   
%   space variable SYSIN to the time domain defined by TIMESTEP. The  
%   optional string argument SYSNAME is used for on-screen messages.

if sysIn.Ts == timeStep
    sysOut = sysIn;
else
    if nargin < 4
        sysName = 'System';
    end
    if sysIn.Ts > 0 && timeStep == 0
        sysOut = d2c(sysIn,method);
    elseif sysIn.Ts == 0 && timeStep > 0
        sysOut = c2d(sysIn,timeStep,method);
    %The following case may not be necessary - using d2d gives the same results    
    elseif sysIn.Ts <= 0.5*timeStep
        sysOut = d2c(sysIn,method);
        sysOut = c2d(sysOut,timeStep,method);
    else
        sysOut = d2d(sysIn,timeStep,method);
    end
end