close all
warning off

run ModalInfo

%% Filter File %%

[Modalinfo,Modalname] = xlsread('.\Temp\Modal.xlsx','Sheet3');

isv47 = 0;
for i = 1:length(Modalname(:,1))
    if size(Modalname{i,1},2) >= 5
        if strcmp(Modalname{i,1}(end-5:end),'cyclic')
            isv47 = 1;
            break
        end
    end
end

Modalname(1,:) = [];%%

if isv47
    %%%%%%%%
    pitch_filter(1,:) = [0,0,10,0.6];
    Pitch_filterlist =  {'Rotor 1st edgewise collective','Rotor 2nd edgewise collective','Tower 2nd side-side mode','Rotor 3rd edgewise collective','3P','6P','9P'};
    pitch_d1 = ones(length(Pitch_filterlist),1);
    pitch_d2 = ones(length(Pitch_filterlist),1);
    %%%%%%%%
    for i = 1:length(Pitch_filterlist)
        pitch_d1(i) = 0.02;
        pitch_d2(i) = 0.5;
    end

    k = 2;
    for i = 1:length(Modalname)
        if ismember(Modalname{i,1}, Pitch_filterlist)
            pitch_filter(k, 1) = Modalinfo(i,1)*2*pi;
            pitch_filter(k, 2) = pitch_d1(k-1);
            pitch_filter(k, 3) = Modalinfo(i,1)*2*pi;
            pitch_filter(k, 4) = pitch_d2(k-1);
            k = k + 1;

            if k > length(Pitch_filterlist) + 1
                break
            end
        end
    end
    %%%%%%%
    torque_filter(1,:) = [0,0,7.87,1];
    Torque_filterlist = {'Tower 1st fore-aft mode','Rotor 1st edgewise sine cyclic','Tower 3rd side-side mode','Rotor 2nd edgewise cosine cyclic','3P','6P'};
    torque_d1 = ones(length(Torque_filterlist),1);
    torque_d2 = ones(length(Torque_filterlist),1);
    %%%%%%%
    for i = 1:length(Torque_filterlist)
        torque_d1(i) = 0.02;
        torque_d2(i) = 0.5;
    end

    k = 2;
    for i = 1:length(Modalname)
        if ismember(Modalname{i,1}, Torque_filterlist)
            torque_filter(k, 1) = Modalinfo(i,1)*2*pi;
            torque_filter(k, 2) = torque_d1(k-1);
            torque_filter(k, 3) = Modalinfo(i,1)*2*pi;
            torque_filter(k, 4) = torque_d2(k-1);
            k = k + 1;

            if k > length(Torque_filterlist) + 1
                break
            end
        end
    end
    %%%%%%
    naf_filter(1,:) = [0,0,10,0.6];
    NAF_filterlist = {'Rotor speed (1P)','3P','6P'};
    naf_d1 = ones(length(NAF_filterlist),1);
    naf_d2 = ones(length(NAF_filterlist),1);
    %%%%%%
    for i = 1:length(NAF_filterlist)
        naf_d1(i) = 0.02;
        naf_d2(i) = 0.5;
    end

    k = 2;
    for i = 1:length(Modalname)
        if ismember(Modalname{i,1}, NAF_filterlist)
            naf_filter(k, 1) = Modalinfo(i,1)*2*pi;
            naf_filter(k, 2) = naf_d1(k-1);
            naf_filter(k, 3) = Modalinfo(i,1)*2*pi;
            naf_filter(k, 4) = naf_d2(k-1);
            k = k + 1;

            if k > length(NAF_filterlist) + 1
                break
            end
        end
    end
else
    %%%%%%
    pitch_filter(1,:) = [0,0,10,0.6];
    Pitch_filterlist =  {'Tower mode 1','Rotor mode 2  B','3P','6P'};
    pitch_d1 = ones(length(Pitch_filterlist),1);
    pitch_d2 = ones(length(Pitch_filterlist),1);
    %%%%%%
    for i = 1:length(Pitch_filterlist)
        pitch_d1(i) = 0.02;
        pitch_d2(i) = 0.5;
    end

    k = 2;
    for i = 1:length(Modalname)
        if ismember(Modalname{i,1}, Pitch_filterlist)
            pitch_filter(k, 1) = Modalinfo(i,1);
            pitch_filter(k, 2) = pitch_d1(k-1);
            pitch_filter(k, 3) = Modalinfo(i,1);
            pitch_filter(k, 4) = pitch_d2(k-1);
            k = k + 1;

            if k > length(Pitch_filterlist) + 1
                break
            end
        end
    end
    %%%%%%%
    torque_filter(1,:) = [0,0,7.87,1];
    Torque_filterlist = {'Tower mode 1','Rotor mode 2  C','3P','6P'};
    torque_d1 = ones(length(Torque_filterlist),1);
    torque_d2 = ones(length(Torque_filterlist),1);
    %%%%%%%
    for i = 1:length(Torque_filterlist)
        torque_d1(i) = 0.02;
        torque_d2(i) = 0.5;
    end

    k = 2;
    for i = 1:length(Modalname)
        if ismember(Modalname{i,1}, Torque_filterlist)
            torque_filter(k, 1) = Modalinfo(i,1);
            torque_filter(k, 2) = torque_d1(k-1);
            torque_filter(k, 3) = Modalinfo(i,1);
            torque_filter(k, 4) = torque_d2(k-1);
            k = k + 1;

            if k > length(Torque_filterlist) + 1
                break
            end
        end
    end
    %%%%%%
    naf_filter(1,:) = [0,0,10,0.6];
    NAF_filterlist = {'1P','3P','6P'};
    naf_d1 = ones(length(NAF_filterlist),1);
    naf_d2 = ones(length(NAF_filterlist),1);
    %%%%%%
    for i = 1:length(NAF_filterlist)
        naf_d1(i) = 0.02;
        naf_d2(i) = 0.5;
    end

    k = 2;
    for i = 1:length(Modalname)
        if ismember(Modalname{i,1}, NAF_filterlist)
            naf_filter(k, 1) = Modalinfo(i,1);
            naf_filter(k, 2) = naf_d1(k-1);
            naf_filter(k, 3) = Modalinfo(i,1);
            naf_filter(k, 4) = naf_d2(k-1);
            k = k + 1;

            if k > length(NAF_filterlist) + 1
                break
            end
        end
    end
end

if exist('.\Temp\Filters.xlsx','file')
    delete('.\Temp\Filters.xlsx');
end

status_1 = xlswrite('.\Temp\Filters.xlsx',pitch_filter,'Pitch','A1');
status_2 = xlswrite('.\Temp\Filters.xlsx',torque_filter,'Torque','A1');
status_3 = xlswrite('.\Temp\Filters.xlsx',naf_filter,'NAF','A1');

if ~status_1 || ~status_2 || ~status_3
    error('Filter.xlsx saved failed')
end

%% Stations %%

load('.\LinearModel\linmod1.mat');

Order = 1;  %Station Number
S_sign = 1;
FinePitch = PitchAngles(1); %Minimum Fine PitchAngles
PitchLoop = 0; %PitchAngles Loop Control Flag


while ~PitchLoop&&S_sign < length(PitchAngles)
    if PitchAngles(S_sign) > FinePitch
        Q_station = S_sign - 1; %% Torque loop
        Q_pitchangle = PitchAngles(S_sign - 1);
        Q_windspeed = Windspeeds(S_sign - 1);
        
        P_station(Order) = S_sign; %%Pitch loop first station
        P_pitchangle(Order) = PitchAngles(S_sign);
        P_windspeed(Order) = Windspeeds(S_sign);
        
        PitchLoop = 1;
    end
    S_sign = S_sign + 1;
end

Order = Order + 1;
P_station(Order) = S_sign + 1;
P_pitchangle(Order) = PitchAngles(S_sign + 1);  %%Second Station
P_windspeed(Order) = Windspeeds(S_sign + 1);


for j = S_sign + 3:4:length(PitchAngles)
    Order = Order + 1;
    if j < (length(PitchAngles)-4)
        P_station(Order) = j;
        P_pitchangle(Order) = PitchAngles(j);  %%Further Station
        P_windspeed(Order) = Windspeeds(j);
    else
        P_station(Order) = length(PitchAngles);
        P_pitchangle(Order) = PitchAngles(end);  %%Further Station
        P_windspeed(Order) = Windspeeds(end);
    end
end

if exist('.\Temp\Station.xlsx','file')
    delete('.\Temp\Station.xlsx');
end

status_1 = xlswrite('.\Temp\Station.xlsx',P_station','Pitch','A1');
status_2 = xlswrite('.\Temp\Station.xlsx',P_pitchangle','Pitch','B1');
status_3 = xlswrite('.\Temp\Station.xlsx',P_windspeed','Pitch','C1');

status_4 = xlswrite('.\Temp\Station.xlsx',Q_station,'Torque','A1');
status_5 = xlswrite('.\Temp\Station.xlsx',Q_pitchangle,'Torque','B1');
status_6 = xlswrite('.\Temp\Station.xlsx',Q_windspeed,'Torque','C1');

if ~status_1 || ~status_2 || ~status_3 || ~status_4 || ~status_5 || ~status_6
    error('Station.xlsx saved failed')
end
