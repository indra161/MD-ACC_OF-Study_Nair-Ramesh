load('B2D1T1_Obs_Freezing_Updated.mat','freezingResponses')

totalHabituation = 0;
totalOFF = 0;
totalON = 0;

freezeLength = length(freezingResponses(:,2));
for i = 1:freezeLength
    if freezingResponses{i,2} < 300
        totalHabituation = totalHabituation + freezingResponses{i,4};
    elseif freezingResponses{i,2} >= 300 && freezingResponses{i,2} < 730
        totalOFF = totalOFF + freezingResponses{i,4};
    else
        totalON = totalON + freezingResponses{i,4};
    end
end

totalHabituation
totalOFF
totalON

% load('B2D1T1_Obs_Freezing_Updated.mat','offFreezeResultsTable', 'onFreezeResultsTable')
% 
% offFreq = 0;
% counter = 0;
% 
% for i = 1:10
%     cellVal = string(offFreezeResultsTable.FreezingOccurred(i));
%     if cellVal == 'Yes'
%         counter = counter + 1;
%     end 
% end
% 
% offFreq = (counter/10)*100;
% 
% onFreq = 0;
% counter = 0;
% 
% for i = 1:10
%     cellVal = string(onFreezeResultsTable.FreezingOccurred(i));
%     if cellVal == 'Yes'
%         counter = counter + 1;
%     end 
% end
% 
% onFreq = (counter/10)*100;
