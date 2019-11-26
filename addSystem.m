function sysOut = addSystem(sys1,sys2,connectionType,sys1Name,sys2Name)
%Adds an additional input/output system in series to a plant.
%
%   SYSOUT = ADDSYSTEM(SYS1,SYS2,CONNECTIONTYPE,SYS1NAME,SYS2NAME) joins   
%   the dynamic systems SYS1 and SYS2 into one system SYSOUT. The input    
%   string CONNECTIONTYPE specifies whether the two systems are connected   
%   in series or in feedback. The function checks the consistency of all  
%   the input and output names of SYS1 and SYS2 and joins them using the
%   index-based approach of the CONNECT function. The external inputs of 
%   SYSOUT are the inputs of SYS1 plus the inputs of SYS2 that do not match 
%   any input or output of SYS1 (in the latter case they are connected to  
%   that output and suppressed). The external outputs of SYSOUT are the  
%   outputs of SYS1 plus the outputs of SYS2, except those that match any 
%   input of SYS1 (this is only allowed in feedback connections). The
%   optional inputs SYS1NAME and SYS2NAME are used for on-screen messages.

if (~strcmp(connectionType,'series') && ~strcmp(connectionType,'feedback'))
    error('connectionType must be either ''series'' or ''feedback''');
end

if nargin < 4
    sys1Name = 'System 1';
end

if nargin < 5
    sys2Name = 'System 2';
end
%==========================================================================
%Check that sys1 and sys2 have the same timesteps
if sys1.Ts ~= sys2.Ts
    error([sys1Name ' and ' sys2Name ' have different timesteps']);
end
%==========================================================================
if (size(sys1.b,2) == 1 && size(sys1.c,1) == 1 && size(sys2.b,2) == 1 && size(sys2.c,1) == 1)
    if strcmp(connectionType,'series')
        sysOut = series(sys1,sys2);
    else
        sysOut = feedback(sys1,sys2);
    end
%==========================================================================
else
    %Check consistency of sys1 and sys2 inputs and outputs for
    %name-based use of connect function
    
    %Check that all sys1 inputs are named
    if sum(strcmp('',cellstr(char(sys1.inputname)))) > 0
        error(['The input names of ' sys1Name ' have to be fully defined.']);
    end
    %Check that all sys1 outputs are named
    if sum(strcmp('',cellstr(char(sys1.outputname)))) > 0
        error(['The output names of ' sys1Name ' have to be fully defined.']);
    end
    %Check that all sys2 inputs are named
    if sum(strcmp('',cellstr(char(sys2.inputname)))) > 0
        error(['The input names of ' sys2Name ' have to be fully defined.']);
    end
    %Check that all sys2 outputs are named
    if sum(strcmp('',cellstr(char(sys2.outputname)))) > 0
        error(['The output names of ' sys2Name ' have to be fully defined.'])
    end

    %For each sys1 input...
    for i = 1:length(sys1.inputname)
        %...check that no other sys1 input has the same name
        if sum(strcmp(sys1.inputname(i),sys1.inputname)) > 1
            error([sys1Name ' has more than one input named ' char(sys1.inputname(i)) '!']);
        end
        %...check that no sys1 output has the same name
        if sum(strcmp(sys1.inputname(i),sys1.outputname)) > 0
            error([sys1Name ' has both an input and an output named ' char(sys1.inputname(i)) '!']);
        end
    end
    
    %For each sys1 output...
    for i = 1:length(sys1.outputname)
        %...check that no other sys1 output has the same name
        if sum(strcmp(sys1.outputname(i),sys1.outputname)) > 1
            error([sys1Name ' has more than one output named ' char(sys1.outputname(i)) '!'])
        end
        %...check that no sys2 output has the same name
        if sum(strcmp(sys1.outputname(i),sys2.outputname)) > 0
            error([sys1Name ' and ' sys2Name ' both have an output named ' char(sys2.outputname(i)) '!'])
        end
    end     
    
    %For each sys2 input...
    for i = 1:length(sys2.inputname)
        %...check that no other sys2 input has the same name
        if sum(strcmp(sys2.inputname(i),sys2.inputname)) > 1
            error([sys2Name ' has more than one input named ' char(sys2.inputname(i)) '!']);
        end
        %...check that no sys2 output has the same name
        if sum(strcmp(sys2.inputname(i),sys2.outputname)) > 0
            error([sys2Name ' has both an input and an output named ' char(sys2.inputname(i)) '!']);
        end
    end    
    
    %For each sys2 output...
    for i = 1:length(sys2.outputname)
        %...check that no other sys2 output has the same name
        if sum(strcmp(sys2.outputname(i),sys2.outputname)) > 1
            error([sysName ' has more than one output named ' char(sys2.outputname(i)) '!'])
        end
        %If the connection is in series... 
        if strcmp(connectionType,'series')
            %...also check that no sys1 input has the same name
            if sum(strcmp(sys2.outputname(i),sys1.inputname)) > 0
                error('For series connections no input of system 1 can match an output of system 2!')
            end
        end        
    end
    %======================================================================    
    %Join sys1 and sys1 in parallel (to combine common inputs)
    sysOut = parallel(sys1,sys2,'names');
    
%     %Combine common inputs
%     appendedInputname = sysOut.inputname;
%     appendedOutputname = sysOut.outputname;
%     b = sysOut.b;
%     d = sysOut.d;
%     i = 1;
%      
%     while i < length(appendedInputname)
%         j = i + 1;
%         while j <= length(appendedInputname)
%             if strcmp(appendedInputname(i),appendedInputname(j))
%                 b(:,i) = b(:,i) + b(:,j);
%                 d(:,i) = d(:,i) + d(:,j);
%                 if j < length(appendedInputname)
%                     b = [b(:,1:j-1) b(:,j+1:end)];                
%                     d = [d(:,1:j-1) d(:,j+1:end)];
%                     appendedInputname = [appendedInputname(1:j-1); appendedInputname(j+1:end)];
%                 else
%                     b = b(:,1:j-1);                
%                     d = d(:,1:j-1);                   
%                     appendedInputname = appendedInputname(1:j-1);        
%                 end
%             else
%                 j = j + 1;
%             end                           
%         end
%         i = i + 1;
%     end
%     
%     sysOut = ss(sysOut.a,b,sysOut.c,d);
%     sysOut.inputname = appendedInputname;
%     sysOut.outputname = appendedOutputname;
    
    %======================================================================      
    %Connect sys1 outputs to sys2 inputs (series connection)
    Q = [];
    newInputIndex = 1:length(sysOut.inputname);
    newOutputIndex = 1:length(sysOut.outputname);
    
    for i = 1:length(sys1.outputname)
        j = find(strcmp(sys1.outputname(i),sysOut.inputname));
        if ~isempty(j)
            newInputIndex(j) = 0;
            Q = [Q; j i];
        end
    end
    
    newInputIndex = newInputIndex(newInputIndex>0);
       
    if isempty(Q)
        error(['At least one output of ' sys1Name ' must match an input of ' sys2Name '!']);
    end
    
    sysOut = connect(sysOut,Q,newInputIndex,newOutputIndex);
    
    %======================================================================     
    %Connect sys2 outputs to sys1 inputs (feedback connection)
    if strcmp(connectionType,'feedback')
        Q = [];
        newInputIndex = 1:length(sysOut.inputname);
        newOutputIndex = 1:length(sysOut.outputname);
        
        for i = 1:length(sys1.inputname)
            j = find(strcmp(sys1.inputname(i),sysOut.outputname));
            if ~isempty(j)
                newOutputIndex(j) = 0;
                Q = [Q; i j];
            end
        end
        
        newOutputIndex = newOutputIndex(newOutputIndex>0);
        
        if isempty(Q)
            error(['For feedback connections at least one output of ' sys2Name ' must match an input of ' sys1Name '!']);
        end  
        
        sysOut = connect(sysOut,Q,newInputIndex,newOutputIndex);
    end     
end
%==========================================================================
