function sysOut = addController(sysIn,controller,ControllerType)
%Adds a controller to a plant in positive feedback.
%   
%   SYSOUT = ADDCONTROLLER(SYSIN,CONTROLLER,CONTROLLERTYPE) returns the 
%   state space model obtained by adding CONTROLLER in positive feedback
%   to SYSIN. 

%Convert plant to the controller type timestep using zero order hold
sysIn = convertSsTime(sysIn,ControllerType.TimeStep,'ZOH','Plant');
%Convert controller to the controller type timestep using Tustin transform
controller = convertSsTime(controller,ControllerType.TimeStep,'tustin',ControllerType.Name);

%Add input and output delays to controller
if ControllerType.TimeStep
    controller.inputdelay = round(ControllerType.InputDelay/ControllerType.TimeStep);
    %Make sure that roundings of input and output delay don't exceed unit
    controller.outputdelay = round((ControllerType.OutputDelay+ControllerType.InputDelay)/ControllerType.TimeStep) - controller.inputdelay;
end

%Add controller to plant
sysOut = addSystem(sysIn,controller,'feedback','Plant',ControllerType.Name);