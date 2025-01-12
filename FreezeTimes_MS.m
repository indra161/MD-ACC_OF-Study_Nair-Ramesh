load('B2D1T1_Obs_Freezing_Updated.mat','freezingResponses')

%% Preshock Freezing

OFF_shock=[300 350 390 420 480 530 560 600 660 690];  
ON_shock=[730 780 820 880 910 950 990 1050 1100 1160];

OFF_PreShock = OFF_shock - 5;
ON_PreShock = ON_shock - 5;

freezeLength = length(freezingResponses(:,2));

%Off Shocks
for i = 1:10
    j = 0;
    for j = 1:freezeLength
        %duration = 0;
        if freezingResponses{j,2} <= OFF_PreShock(i) && OFF_PreShock(i) <= freezingResponses{j,3} %%freezing starts before preshock time
            if freezingResponses{j,3} >= OFF_shock(i)
                freezingResponses{j,5} = OFF_shock(i) - OFF_PreShock(i); %%if freezing ends during shock period, preshock freeze = 5 sec
            else
                freezingResponses{j,5} = freezingResponses{j,3} - OFF_PreShock(i); %%if freezing ends before shock period, preshock freeze = end time - preshcok start
            end
        elseif freezingResponses{j,2} >= OFF_PreShock(i) && freezingResponses{j,2} < OFF_shock(i) %freezing starts after preshock period
            if freezingResponses{j,3} >= OFF_shock(i)
                freezingResponses{j,5} = OFF_shock(i) - freezingResponses{j,2}; %%if freezing ends during shock period, preshock freeze = shock time - preshock start
            else
                freezingResponses{j,5} = freezingResponses{j,3} - freezingResponses{j,2}; %%if freezing ends before shock period, preshock freeze = end time - preshcok start
            end
        % else
        %     duration = 0;
        end
        %freezingResponses{j,5} = duration;
    end
end 

%ON Shocks
for i = 1:10
    j = 0;
    for j = 1:freezeLength
        %duration = 0;
        if freezingResponses{j,2} <= ON_PreShock(i) && ON_PreShock(i) <= freezingResponses{j,3} %%freezing starts before preshock time
            if freezingResponses{j,3} >= ON_shock(i)
                freezingResponses{j,5} = ON_shock(i) - ON_PreShock(i); %%if freezing ends during shock period, preshock freeze = 5 sec
            else
                freezingResponses{j,5} = freezingResponses{j,3} - ON_PreShock(i); %%if freezing ends before shock period, preshock freeze = end time - preshcok start
            end
        elseif freezingResponses{j,2} >= ON_PreShock(i) && freezingResponses{j,2} < ON_shock(i) %freezing starts after preshock period
            if freezingResponses{j,3} >= ON_shock(i)
                freezingResponses{j,5} = ON_shock(i) - freezingResponses{j,2}; %%if freezing ends during shock period, preshock freeze = shock time - preshock start
            else
                freezingResponses{j,5} = freezingResponses{j,3} - freezingResponses{j,2}; %%if freezing ends before shock period, preshock freeze = end time - preshcok start
            end
        % else
        %     duration = 0;
        end
        %freezingResponses{j,5} = duration;
    end
end 

%convert empty cells to zero
for i = 1:freezeLength
    if isempty(freezingResponses{i,5})
        freezingResponses{i,5} = 0;
    end
end

%% PostShock Period Freezing (Shock Period + 3 seconds post shock --> total 5s from start of shock)

OFF_PostShock = OFF_shock + 5;
ON_PostShock = ON_shock + 5;

%Off Shocks
for i = 1:10
    j = 0;
    for j = 1:freezeLength
        if freezingResponses{j,2} >= OFF_shock(i) && freezingResponses{j,2} <= OFF_PostShock(i) %%freezing starts between shock moment and end of postshock moment
            if freezingResponses{j,3} >= OFF_PostShock(i)
                freezingResponses{j,6} = OFF_PostShock(i) - freezingResponses{j,2}; %%if freezing ends after post shock period, time = postshock end - start time
            else
                freezingResponses{j,6} = freezingResponses{j,3} - freezingResponses{j,2}; %%if freezing end before postshock end, time = end - start
            end
        elseif freezingResponses{j,2} <= OFF_shock(i) && freezingResponses{j,3} >= OFF_shock(i) %%if freezing starts before shock moment but continues into postshock period, time = end time - shock moment
            if freezingResponses{j,3} >= OFF_PostShock(i)
                freezingResponses{j,6} = OFF_PostShock(i) - OFF_shock(i); %%if freezing ends after post shock period, time = postshock end - start time
            else
                freezingResponses{j,6} = freezingResponses{j,3} - OFF_shock(i); %%if freezing end before postshock end, time = end - start
            end
        end
    end
end 

%ON Shocks
for i = 1:10
    j = 0;
    for j = 1:freezeLength
        if freezingResponses{j,2} >= ON_shock(i) && freezingResponses{j,2} <= ON_PostShock(i) %%freezing starts between shock moment and end of postshock moment
            if freezingResponses{j,3} >= ON_PostShock(i)
                freezingResponses{j,6} = ON_PostShock(i) - freezingResponses{j,2}; %%if freezing ends after post shock period, time = postshock end - start time
            else
                freezingResponses{j,6} = freezingResponses{j,3} - freezingResponses{j,2}; %%if freezing end before postshock end, time = end - start
            end
        elseif freezingResponses{j,2} <= ON_shock(i) && freezingResponses{j,3} >= ON_shock(i) %%if freezing starts before shock moment but continues into postshock period, time = end time - shock moment
            if freezingResponses{j,3} >= ON_PostShock(i)
                freezingResponses{j,6} = ON_PostShock(i) - ON_shock(i); %%if freezing ends after post shock period, time = postshock end - start time
            else
                freezingResponses{j,6} = freezingResponses{j,3} - ON_shock(i); %%if freezing end before postshock end, time = end - start
            end
        end
    end
end 

%convert empty cells to zero
for i = 1:freezeLength
    if isempty(freezingResponses{i,6})
        freezingResponses{i,6} = 0;
    end
end

%% Shock Freezing Times (2 sec shock moment)

OFF_ShockEnd = OFF_shock + 2;
ON_ShockEnd = ON_shock +2;

%Off Shocks
for i = 1:10
    j = 0;
    for j = 1:freezeLength
        if freezingResponses{j,2} >= OFF_shock(i) && freezingResponses{j,2} <= OFF_ShockEnd(i) %%freezing starts between shock moment and end of shock moment
            if freezingResponses{j,3} >= OFF_ShockEnd(i)
                freezingResponses{j,7} = OFF_ShockEnd(i) - freezingResponses{j,2}; %%if freezing ends after shock period, time = shock end - start time
            else
                freezingResponses{j,7} = freezingResponses{j,3} - freezingResponses{j,2}; %%if freezing end before shock end, time = end - start
            end
        elseif freezingResponses{j,2} <= OFF_shock(i) && freezingResponses{j,3} >= OFF_shock(i) %%if freezing starts before shock moment but continues into shock period
            if freezingResponses{j,3} >= OFF_ShockEnd(i)
                freezingResponses{j,7} = OFF_ShockEnd(i) - OFF_shock(i); %%if freezing ends after shock period, time = shock end - shock start
            else
                freezingResponses{j,7} = freezingResponses{j,3} - OFF_shock(i); %%if freezing end before shock end, time = end - shock start
            end
        end
    end
end 

%On Shocks
for i = 1:10
    j = 0;
    for j = 1:freezeLength
        if freezingResponses{j,2} >= ON_shock(i) && freezingResponses{j,2} <= ON_ShockEnd(i) %%freezing starts between shock moment and end of shock moment
            if freezingResponses{j,3} >= ON_ShockEnd(i)
                freezingResponses{j,7} = ON_ShockEnd(i) - freezingResponses{j,2}; %%if freezing ends after shock period, time = shock end - start time
            else
                freezingResponses{j,7} = freezingResponses{j,3} - freezingResponses{j,2}; %%if freezing end before shock end, time = end - start
            end
        elseif freezingResponses{j,2} <= ON_shock(i) && freezingResponses{j,3} >= ON_shock(i) %%if freezing starts before shock moment but continues into shock period
            if freezingResponses{j,3} >= ON_ShockEnd(i)
                freezingResponses{j,7} = ON_ShockEnd(i) - ON_shock(i); %%if freezing ends after shock period, time = shock end - shock start
            else
                freezingResponses{j,7} = freezingResponses{j,3} - ON_shock(i); %%if freezing end before shock end, time = end - shock start
            end
        end
    end
end 

%convert empty cells to zero
for i = 1:freezeLength
    if isempty(freezingResponses{i,7})
        freezingResponses{i,7} = 0;
    end
end
%% Total Freezing Times


totalPreFreeze_OFF = 0;


for i = 1:freezeLength
    if freezingResponses{i,2} < ON_PreShock(1)
        totalPreFreeze_OFF = totalPreFreeze_OFF + freezingResponses{i,5};
    end
end

totalPostFreeze_OFF = 0;

for i = 1:freezeLength
    if freezingResponses{i,2} < ON_PostShock(1)
    totalPostFreeze_OFF = totalPostFreeze_OFF + freezingResponses{i,6};
    end
end

totalShockFreeze_OFF = 0;
for i = 1:freezeLength
    if freezingResponses{i,2} < ON_shock(1)
    totalShockFreeze_OFF = totalShockFreeze_OFF + freezingResponses{i,7};
    end
end

totalPreFreezePercent_OFF = (totalPreFreeze_OFF/50)*100;
totalPostFreezePercent_OFF = (totalPostFreeze_OFF/50)*100;
totalShockFreezePercent_OFF = (totalShockFreeze_OFF/20)*100;

disp('Total freeze time in OFF Preshock Period:')
disp(totalPreFreeze_OFF)
disp('Total freeze time in OFF Postshock Period:')
disp(totalPostFreeze_OFF)
disp('Total freeze time in OFF Shock Period:')
disp(totalShockFreeze_OFF)
disp('Total freeze time in OFF Preshock Period (%):')
disp(totalPreFreezePercent_OFF)
disp('Total freeze time in OFF Postshock Period (%):')
disp(totalPostFreezePercent_OFF)
disp('Total freeze time in OFF Shock Period (%):')
disp(totalShockFreezePercent_OFF)


totalPreFreeze_ON = 0;


for i = 1:freezeLength
    if freezingResponses{i,2} > OFF_PreShock(10)
        totalPreFreeze_ON = totalPreFreeze_ON + freezingResponses{i,5};
    end
end

totalPostFreeze_ON = 0;

for i = 1:freezeLength
    if freezingResponses{i,2} > OFF_PostShock(10)
    totalPostFreeze_ON = totalPostFreeze_ON + freezingResponses{i,6};
    end
end

totalShockFreeze_ON = 0;

for i = 1:freezeLength
    if freezingResponses{i,2} > OFF_shock(10)
    totalShockFreeze_ON = totalShockFreeze_ON + freezingResponses{i,7};
    end
end

totalPreFreezePercent_ON = (totalPreFreeze_ON/50)*100;
totalPostFreezePercent_ON = (totalPostFreeze_ON/50)*100;
totalShockFreezePercent_ON = (totalShockFreeze_ON/20)*100;

disp('Total freeze time in ON Preshock Period:')
disp(totalPreFreeze_ON)
disp('Total freeze time in ON Postshock Period:')
disp(totalPostFreeze_ON)
disp('Total freeze time in ON Shock Period:')
disp(totalShockFreeze_ON)
disp('Total freeze time in ON Preshock Period (%):')
disp(totalPreFreezePercent_ON)
disp('Total freeze time in ON Postshock Period (%):')
disp(totalPostFreezePercent_ON)
disp('Total freeze time in ON Shock Period (%):')
disp(totalShockFreezePercent_ON)

%% Avg Number of freezing epochs

numOFF = 0;

for i = 1:freezeLength
    if freezingResponses{i,2} < ON_shock(1) && freezingResponses{i,2} >= OFF_shock(1)
       numOFF = numOFF + 1;
    end
end

numON = 0;

for i = 1:freezeLength
    if freezingResponses{i,2} > OFF_shock(10)
        numON = numON + 1;
    end
end
