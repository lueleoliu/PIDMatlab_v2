function sysOut = getSiso(sysIn,input,output)
%Extracts a SISO state space model from a MIMO system.
%   
%   SYSOUT = GETSISO(SYSIN,INPUT,OUTPUT) extracts from the MIMO system 
%   SYSIN the SISO state space model from INPUT to OUTPUT.

iInput = find((strcmp(input,sysIn.inputname)));
iOutput = find((strcmp(output,sysIn.outputname)));

%Check that at most one system input has the desired name
if length(iInput) > 1
    error(['System has more than one input named ' input '!'])
end
%Check that at least one system input has the desired name
if length(iInput) < 1
    error(['System has no inputs named ' input '!'])
end
%Check that at most one system output has the desired name
if length(iOutput) > 1
    error(['System has more than one output named ' output '!'])
end
%Check that at least one system output has the desired name
if length(iOutput) < 1
    error(['System has no outputs named ' output '!'])
end        
        
sysOut = ss(sysIn(iOutput,iInput));
sysOut.inputname = input;
sysOut.outputname = output;