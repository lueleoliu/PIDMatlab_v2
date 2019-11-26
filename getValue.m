function value = getValue(name,config)
%This function is defined to return the value of specific config parameter
msgID_2 = 'SettingError:';
% Check config format %
msgID_2_1 = [msgID_2 , 'FormatError'];

if ~iscell(config)
    errmsg_2_1_1 = 'Error:Error occurs during setting uploading, please check the setting file and upload again';
    error(msgID_2_1,errmsg_2_1_1);
end

if size(config,2) > 2
    errmsg_2_1_2 = 'Error:There is unexpected info in config page, only two columns are allowed, please please check the setting file and upload again ';
    error(msgID_2_1,errmsg_2_1_2);
end

set_name = config(:,1);
try
    set_value = cell2mat(config(:,2));
catch ME
    if (strcmp(ME.identifier,'MATLAB:cell2mat:MixedDataTypes'))
        errmsg_2_1_3 = 'Error:Illegal value in the config page, the second column should be all in numbers, please check the setting file and upload again';
        error(msgID_2_1,errmsg_2_1_3);
    else
        rethrow(ME)
    end
end

if size(set_value) ~= size(set_name)
    errmsg_2_1_4 = 'Error:The parameters do not match there values in config page, please check the setting file and upload again ';
    error(msgID_2_1,errmsg_2_1_4);
end

[exist, pos] = ismember(name, set_name);

if ~exist
    errmsg_2_1_5 = ['Error:',name, ' is required to be set(This parameter is not listed), please check the setting file and upload again '];
    error(msgID_2_1,errmsg_2_1_5);
end

r_value = set_value(pos);

if isnan(r_value)
    errmsg_2_1_6 = ['Error:',name, ' value is missing, please check the setting file and upload again '];
    error(msgID_2_1,errmsg_2_1_6);
end

value = r_value;

end

