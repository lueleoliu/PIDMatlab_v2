close all
warning off

%% Turbine Parameters %%
%Default Parameters%
file_i = fopen('.\Campbell\Campbell.%02');
if file_i == -1
    error('File read failed')
end

while ~feof(file_i)
    str = fgetl(file_i);
    if ~isempty(str)
        a = regexp(str,'\s');
        
        if strcmp(str(1:a(1)-1), 'NDIMENS')
            nd = str2num(str(a:end));
        end
        
        if strcmp(str(1:a(1)-1), 'ACCESS')
            acc = str(a+1:end);
        end
        
        if strcmp(str(1:a(1)-1), 'AXITICK')
            a = regexp(str,'''');
            j = 1;
            for i = 1:2:length(a)
                Modal{j} = str(a(i)+1:a(i+1)-1);
                j = j+1;
            end
            if ~strcmp(acc,'D') || nd ~= 3
                Modal_R = {'1P','2P','3P','4P','5P','6P','7P','8P','9P'};
                Modal = [Modal,Modal_R];
            end
        end
      
        if strcmp(str(1:a(1)-1), 'AXIVAL')
            z = 2;
            if strcmp(acc,'D') && nd == 3
                z = 1;
            end
            for i = z:length(a)-1
                Speed{i-z+1} = str2num(str(a(i)+1:a(i+1)-1));
            end
            Speed{i-z+2} = str2num(str(a(i+1)+1:end));
        end
    end
end
fclose(file_i);


%Inherent Frequency%
if strcmp(acc,'D') && nd == 3
    fid = fopen('.\Campbell\Campbell.$02');
    f = fread(fid,'double');
    j = 1;
    for i = 1:3:length(f)
        A(j,1) = f(i);
        A(j,2) = f(i+1);
        A(j,3) = f(i+2);
        j = j + 1;
    end
    
    Fre_a = A(:,1);
    Damping_a = A(:,2);
    
    position = 1;
    for i = 1:length(Speed)
        for j = 1:length(Modal)
        Fre_f(j,i) = Fre_a(position);
        Damping_f(j,i) = Damping_a(position);
        position = position+1;
        end
    end
    
    Rate_p = length(Speed);
    
else
    
    [Fre_a,Damping_a] = textread('.\Campbell\Campbell.$02','%f %f');

    position = 1;
    for i = 1:length(Speed)
        for j = 1:length(Modal)
        Fre_f(j,i) = Fre_a(position);
        Damping_f(j,i) = Damping_a(position);
        position = position+1;
        end
    end

    Err = 1e-4;
    for i = 1:length(Speed)
        if (abs(Speed{i}-Speed{end}) <= Err) && (abs(Speed{i} - Speed{i+1}) <= Err)
            RatedSpeed = Speed{i};
            Rate_p = i;
            break;
        end
    end

end
        
Fre_r = Fre_f(:,Rate_p);
Damping_r = Damping_f(:,Rate_p);


%% Output for check %%
if exist('.\Temp\Modal.xlsx','file')
    delete('.\Temp\Modal.xlsx');
end


%Frequency sheet%
Label = {'Rotor Speed [rad/s]'};
status_1 = xlswrite('.\Temp\Modal.xlsx',Speed,'Sheet1','B1');
status_2 = xlswrite('.\Temp\Modal.xlsx',Modal','Sheet1','A2');
status_3 = xlswrite('.\Temp\Modal.xlsx',Fre_f,'Sheet1','B2');
status_4 = xlswrite('.\Temp\Modal.xlsx',Label,'Sheet1','A1');

if ~status_1
    fprintf('\nModal File Sheet 1 Rotor Speed Info saved failed!\n');
elseif ~status_2
    fprintf('\nModal File Sheet 1 Modal name Info saved failed!\n');
elseif ~status_3
    fprintf('\nModal File Sheet 1 Modal data saved failed!\n');  
elseif ~status_4
    fprintf('\nModal File Sheet 1 Speed Label saved failed!\n'); 
else
    fprintf('\nModal File Sheet 1 Frequency saved successful!\n');
end

%Damping factor sheet%
status_1 = xlswrite('.\Temp\Modal.xlsx',Speed,'Sheet2','B1');
status_2 = xlswrite('.\Temp\Modal.xlsx',Modal','Sheet2','A2');
status_3 = xlswrite('.\Temp\Modal.xlsx',Damping_f,'Sheet2','B2');
status_4 = xlswrite('.\Temp\Modal.xlsx',Label,'Sheet2','A1');

if ~status_1
    fprintf('\nModal File Sheet 2 Rotor Speed Info saved failed!\n');
elseif ~status_2
    fprintf('\nModal File Sheet 2 Modal name Info saved failed!\n');
elseif ~status_3
    fprintf('\nModal File Sheet 2 Modal data saved failed!\n');  
elseif ~status_4
    fprintf('\nModal File Sheet 2 Speed Label saved failed!\n'); 
else
    fprintf('\nModal File Sheet 2 Damping factor saved successful!\n');
end

%Freq vs. Damping at rated rotor speed%
Label_f = {'Frequency [rad/s]'};
Label_d = {'Damping[-]'};
status_1 = xlswrite('.\Temp\Modal.xlsx',Fre_r,'Sheet3','B2');
status_2 = xlswrite('.\Temp\Modal.xlsx',Damping_r,'Sheet3','C2');
status_3 = xlswrite('.\Temp\Modal.xlsx',Modal','Sheet3','A2');
status_4 = xlswrite('.\Temp\Modal.xlsx',Label_f,'Sheet3','B1');
status_5 = xlswrite('.\Temp\Modal.xlsx',Label_d,'Sheet3','C1');

if ~status_1
    fprintf('\nModal File Sheet 3 Frequency saved failed!\n');
elseif ~status_2
    fprintf('\nModal File Sheet 3 Damping Factor saved failed!\n');
elseif ~status_3
    fprintf('\nModal File Sheet 3 Modal Info saved failed!\n');  
elseif ~status_4||~status_5
    fprintf('\nModal File Sheet 3 Label saved failed!\n'); 
else
    fprintf('\nModal File Sheet 3 Modal at rated saved successful!\n');
end

