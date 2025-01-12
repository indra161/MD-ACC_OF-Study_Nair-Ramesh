

%% 

response_Data = importfile_data_analysis(uigetfile('.csv'),[2, Inf]) % load the freezing CSV Data
regroupedData = regroupBehaviorData(response_Data)
FreezingThreshold= 2 % 2 seconds;
freezingResponses={};

for i = 1:height(regroupedData)
    if strcmp(regroupedData.Behavior{i}, 'Freezing')
        duration = regroupedData.EndTime(i) - regroupedData.StartTime(i);
        
        % Check if the duration is more than 2 seconds
        if duration >FreezingThreshold
            % Append relevant freezing response to the new cell array
            freezingResponses = [freezingResponses; {regroupedData.Behavior{i}, regroupedData.StartTime(i), regroupedData.EndTime(i)}];
        end
    end
end


% Calculate duration for each freezing response and add it to the cell array
for i = 1:size(freezingResponses, 1)
    startTime = freezingResponses{i, 2};   % Freezing start time
    endTime = freezingResponses{i, 3};     % Freezing end time
    duration = endTime - startTime;         % Calculate duration
    freezingResponses{i, 4} = duration;     % Add duration to the cell array
end

load('B2D1T1_Data.mat',['Spk_data']) 

OFF_shock=[300 350 390 420 480 530 560 600 660 690];  
ON_shock=[730 780 820 880 910 950 990 1050 1100 1160];


%% Freezing behavior during Shock Period in both OFF and ON conditions
freezeTimesInOFFShock = {};

% Loop through each freezing response
for i = 1:size(freezingResponses, 1)
    startTime = freezingResponses{i, 2}; % Freezing start time
    
    % Check against OFF_shock periods
    for j = 1:length(OFF_shock)
        if startTime >= OFF_shock(j) && startTime < OFF_shock(j) + 2
            freezeTimesInOFFShock = [freezeTimesInOFFShock; {freezingResponses{i, 1}, startTime, OFF_shock(j)}];
        end
    end
end


% Initialize a cell array to hold results for OFF_shock
offFreezeResults = {};

% Loop through each OFF_shock period
for j = 1:length(OFF_shock)
    offShockStart = OFF_shock(j);               % Start of the OFF_shock period
    offShockEnd = OFF_shock(j) + 2;             % End of the OFF_shock period (2 seconds later)
    offStatus = 'No';                           % Default status
    offDuration = NaN;                          % Default duration is NaN

    % Check against freezing responses
    for i = 1:size(freezingResponses, 1)
        startTime = freezingResponses{i, 2}; % Freezing start time
        endTime = freezingResponses{i, 3};   % Freezing end time
        
        % Check if freezing response overlaps with the OFF_shock period
        if (startTime >= offShockStart && startTime < offShockEnd)
            offStatus = 'Yes';               % Freezing happened in this OFF_shock period
            offDuration = endTime - startTime; % Calculate duration
            break;                            % Exit the loop after finding the first match
        end
    end
    
    % Append results to offFreezeResults cell array
    offFreezeResults = [offFreezeResults; {offShockStart, offShockEnd, offStatus, offDuration}];
end

% Convert offFreezeResults to a table for better readability
offFreezeResultsTable = cell2table(offFreezeResults, 'VariableNames', {'ShockStart', 'ShockEnd', 'FreezingOccurred', 'Duration'});

% Display the results for OFF_shock
disp('Freezing response analysis in OFF_shock periods:');
disp(offFreezeResultsTable);

% Initialize a cell array to hold results for ON_shock
onFreezeResults = {};

% Loop through each ON_shock period
for j = 1:length(ON_shock)
    onShockStart = ON_shock(j);               % Start of the ON_shock period
    onShockEnd = ON_shock(j) + 2;             % End of the ON_shock period (2 seconds later)
    onStatus = 'No';                          % Default status
    onDuration = NaN;                         % Default duration is NaN

    % Check against freezing responses
    for i = 1:size(freezingResponses, 1)
        startTime = freezingResponses{i, 2}; % Freezing start time
        endTime = freezingResponses{i, 3};   % Freezing end time
        
        % Check if freezing response overlaps with the ON_shock period
        if (startTime >= onShockStart && startTime < onShockEnd)
            onStatus = 'Yes';                % Freezing happened in this ON_shock period
            onDuration = endTime - startTime; % Calculate duration
            break;                           % Exit the loop after finding the first match
        end
    end
    
    % Append results to onFreezeResults cell array
    onFreezeResults = [onFreezeResults; {onShockStart, onShockEnd, onStatus, onDuration}];
end

% Convert onFreezeResults to a table for better readability
onFreezeResultsTable = cell2table(onFreezeResults, 'VariableNames', {'ShockStart', 'ShockEnd', 'FreezingOccurred', 'Duration'});

% Display the results for ON_shock
disp('Freezing response analysis in ON_shock periods:');
disp(onFreezeResultsTable);


%% Create the response structure of Freezing in the entire timeline of shock experiment
% Create a figure for the timeline
figure;
hold on;

% Plot OFF_shock periods
for j = 1:length(OFF_shock)
    fill([OFF_shock(j), OFF_shock(j) + 2, OFF_shock(j) + 2, OFF_shock(j)], ...
         [0, 0, 1, 1], 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
end

% Plot ON_shock periods
for j = 1:length(ON_shock)
    fill([ON_shock(j), ON_shock(j) + 2, ON_shock(j) + 2, ON_shock(j)], ...
         [0, 0, 1, 1], 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
end

% Plot freezing responses
for i = 1:size(freezingResponses, 1)
    startTime = freezingResponses{i, 2}; % Freezing start time
    endTime = freezingResponses{i, 3};   % Freezing end time
    plot([startTime, endTime], [0.5, 0.5], 'k', 'LineWidth', 3); % Black line for freezing
end

% Set the axis limits and labels
xlim([0, 1200]);
ylim([-0.1, 1.1]);
xlabel('Time (seconds)');
ylabel('Events');
title('Timeline of OFFshock, ONshock, and Freezing Responses');
legend({'OFFshock', 'ONshock', 'Freezing Responses'}, 'Location', 'NorthEast');
grid on;
hold off;

%% Freezing behavior during PostShock Period in both OFF and ON conditions
% Initialize results for post-shock freezing responses
postFreezeResults_OFF = {};
postFreezeResults_ON = {};

% Calculate post-shock periods and check freezing responses for OFF_shock
for j = 1:length(OFF_shock)
    postShockStart = OFF_shock(j) + 2;             % Start of the post-shock period
    postShockEnd = postShockStart + 5;              % End of the post-shock period (2+5 seconds)
    postStatus = 'No';                             % Default status
    postDuration = 0;                              % Default duration is 0

    % Check against freezing responses
    for i = 1:size(freezingResponses, 1)
        startTime = freezingResponses{i, 2};      % Freezing start time
        endTime = freezingResponses{i, 3};        % Freezing end time
        
        % Check if freezing response overlaps with the post-shock period
        if (startTime >= postShockStart && startTime < postShockEnd)
            postStatus = 'Yes';                   % Freezing happened in this post-shock period
            postDuration = endTime-startTime; % Calculate duration
        end
    end
    
    % Append results to postFreezeResults_OFF cell array
    postFreezeResults_OFF = [postFreezeResults_OFF; {postShockStart, postShockEnd, postStatus, postDuration}];
end

% Calculate post-shock periods and check freezing responses for ON_shock
for j = 1:length(ON_shock)
    postShockStart = ON_shock(j) + 2;              % Start of the post-shock period
    postShockEnd = postShockStart + 5;              % End of the post-shock period (2+5 seconds)
    postStatus = 'No';                             % Default status
    postDuration = 0;                              % Default duration is 0

    % Check against freezing responses
    for i = 1:size(freezingResponses, 1)
        startTime = freezingResponses{i, 2};      % Freezing start time
        endTime = freezingResponses{i, 3};        % Freezing end time
        
        % Check if freezing response overlaps with the post-shock period
        if (startTime >= postShockStart && startTime < postShockEnd)
            postStatus = 'Yes';                   % Freezing happened in this post-shock period
            postDuration = endTime-startTime; % Calculate duration
        end
    end
    
    % Append results to postFreezeResults_ON cell array
    postFreezeResults_ON = [postFreezeResults_ON; {postShockStart, postShockEnd, postStatus, postDuration}];
end

% Convert results to tables for better readability
postFreezeResultsTable_OFF = cell2table(postFreezeResults_OFF, 'VariableNames', {'PostShockStart', 'PostShockEnd', 'FreezingOccurred', 'Duration'});
postFreezeResultsTable_ON = cell2table(postFreezeResults_ON, 'VariableNames', {'PostShockStart', 'PostShockEnd', 'FreezingOccurred', 'Duration'});

% Display the results for post-shock periods
disp('Freezing response analysis in post-shock periods for OFF_shock:');
disp(postFreezeResultsTable_OFF);

disp('Freezing response analysis in post-shock periods for ON_shock:');
disp(postFreezeResultsTable_ON);

%% Freezing behavior during PreShock Period in both OFF and ON conditions
% Initialize results for pre-shock freezing responses
preFreezeResults_OFF = {};
preFreezeResults_ON = {};

% Calculate pre-shock periods and check freezing responses for OFF_shock
for j = 1:length(OFF_shock)
    preShockStart = OFF_shock(j) - 3;            % Start of the pre-shock period
    preShockEnd = OFF_shock(j) - 0.01;            % End of the pre-shock period
    preStatus = 'No';                             % Default status
    preDuration = 0;                              % Default duration is 0

    % Check against freezing responses
    for i = 1:size(freezingResponses, 1)
        startTime = freezingResponses{i, 2};      % Freezing start time
        endTime = freezingResponses{i, 3};        % Freezing end time
        
        % Check if freezing response overlaps with the pre-shock period
        if (startTime < preShockEnd && endTime > preShockStart)
            preStatus = 'Yes';                   % Freezing happened in this pre-shock period
            preDuration = preDuration + (min(endTime, preShockEnd) - max(startTime, preShockStart)); % Calculate duration
        end
    end
    
    % Append results to preFreezeResults_OFF cell array
    preFreezeResults_OFF = [preFreezeResults_OFF; {preShockStart, preShockEnd, preStatus, preDuration}];
end

% Calculate pre-shock periods and check freezing responses for ON_shock
for j = 1:length(ON_shock)
    preShockStart = ON_shock(j) - 3;              % Start of the pre-shock period
    preShockEnd = ON_shock(j) - 0.01;              % End of the pre-shock period
    preStatus = 'No';                             % Default status
    preDuration = 0;                              % Default duration is 0

    % Check against freezing responses
    for i = 1:size(freezingResponses, 1)
        startTime = freezingResponses{i, 2};      % Freezing start time
        endTime = freezingResponses{i, 3};        % Freezing end time
        
        % Check if freezing response overlaps with the pre-shock period
        if (startTime < preShockEnd && endTime > preShockStart)
            preStatus = 'Yes';                   % Freezing happened in this pre-shock period
            preDuration = preDuration + (min(endTime, preShockEnd) - max(startTime, preShockStart)); % Calculate duration
        end
    end
    
    % Append results to preFreezeResults_ON cell array
    preFreezeResults_ON = [preFreezeResults_ON; {preShockStart, preShockEnd, preStatus, preDuration}];
end

% Convert results to tables for better readability
preFreezeResultsTable_OFF = cell2table(preFreezeResults_OFF, 'VariableNames', {'PreShockStart', 'PreShockEnd', 'FreezingOccurred', 'Duration'});
preFreezeResultsTable_ON = cell2table(preFreezeResults_ON, 'VariableNames', {'PreShockStart', 'PreShockEnd', 'FreezingOccurred', 'Duration'});

% Display the results for pre-shock periods
disp('Freezing response analysis in pre-shock periods for OFF_shock:');
disp(preFreezeResultsTable_OFF);

disp('Freezing response analysis in pre-shock periods for ON_shock:');
disp(preFreezeResultsTable_ON);

%%

% Initialize results for categorization
categorizedFreezes_OFF = {};
categorizedFreezes_ON = {};
totalDuration_OFF = struct('preShock', 0, 'shock', 0, 'postShock', 0, 'outOfBound', 0);
totalDuration_ON = struct('preShock', 0, 'shock', 0, 'postShock', 0, 'outOfBound', 0);

% Categorize freezing responses for OFF_shock
for i = 1:size(freezingResponses, 1)
    startTime = freezingResponses{i, 2};   % Freezing start time
    endTime = freezingResponses{i, 3};     % Freezing end time
    isCategorized = false;                 % Flag to check if categorized

    % Check for OFF_shock periods
    for j = 1:length(OFF_shock)
       
        % % Pre-shock period
        % if (startTime < OFF_shock(j) && startTime > (OFF_shock(j) - 3))
        %     categorizedFreezes_OFF = [categorizedFreezes_OFF; {'pre_shock_freeze', startTime, endTime}];
        %     totalDuration_OFF.preShock = totalDuration_OFF.preShock + (endTime-startTime);
        %     isCategorized = true;
        %     break;  % No need to check further
        % end
        
        % Shock period
        if (startTime < OFF_shock(j) + 2 && startTime > OFF_shock(j))
            categorizedFreezes_OFF = [categorizedFreezes_OFF; {'shock_freeze', startTime, endTime}];
            totalDuration_OFF.shock = totalDuration_OFF.shock + (endTime-startTime);
            isCategorized = true;
            break;  % No need to check further
        end
        
        % Post-shock period
        if (startTime < OFF_shock(j) + 7 && startTime > (OFF_shock(j) + 2))
            categorizedFreezes_OFF = [categorizedFreezes_OFF; {'post_shock_freeze', startTime, endTime}];
            totalDuration_OFF.postShock = totalDuration_OFF.postShock + (endTime-startTime);
            isCategorized = true;
            break;  % No need to check further
        end
    end

    % If not categorized in OFF_shock, it's out of bound
    if ~isCategorized
        if startTime <= OFF_shock(end)+30
        categorizedFreezes_OFF = [categorizedFreezes_OFF; {'out_of_bound_freeze', startTime, endTime}];
        totalDuration_OFF.outOfBound = totalDuration_OFF.outOfBound + (endTime - startTime);
        end
    end

end



% Categorize freezing responses for ON_shock
for i = 1:size(freezingResponses, 1)
    startTime = freezingResponses{i, 2};   % Freezing start time
    endTime = freezingResponses{i, 3};     % Freezing end time
    isCategorized = false;                 % Flag to check if categorized

    % Check for ON_shock periods
    for j = 1:length(ON_shock)

        % % Pre-shock period
        % if (startTime < ON_shock(j) && startTime > (ON_shock(j) - 3))
        %     categorizedFreezes_ON = [categorizedFreezes_ON; {'pre_shock_freeze', startTime, endTime}];
        %     totalDuration_ON.preShock = totalDuration_ON.preShock +  (endTime - startTime);
        %     isCategorized = true;
        %     break;  % No need to check further
        % end
        % 
        % Shock period
        if (startTime < ON_shock(j) + 2 && startTime > ON_shock(j))
            categorizedFreezes_ON = [categorizedFreezes_ON; {'shock_freeze', startTime, endTime}];
            totalDuration_ON.shock = totalDuration_ON.shock +  (endTime - startTime);
            isCategorized = true;
            break;  % No need to check further
        end
        
        % Post-shock period
        if (startTime < ON_shock(j) + 7 && startTime > (ON_shock(j) + 2))
            categorizedFreezes_ON = [categorizedFreezes_ON; {'post_shock_freeze', startTime, endTime}];
            totalDuration_ON.postShock = totalDuration_ON.postShock +  (endTime - startTime);
            isCategorized = true;
            break;  % No need to check further
        end
    end

    % If not categorized in ON_shock, it's out of bound
    if ~isCategorized
        if startTime>= ON_shock(1)-5
        categorizedFreezes_ON = [categorizedFreezes_ON; {'out_of_bound_freeze', startTime, endTime}];
        totalDuration_ON.outOfBound = totalDuration_ON.outOfBound + (endTime - startTime);
        end
    end
end

% Convert categorized results for OFF_shock to a table for better readability
categorizedTable_OFF = cell2table(categorizedFreezes_OFF, 'VariableNames', {'FreezeType', 'StartTime', 'EndTime'});

% Convert categorized results for ON_shock to a table for better readability
categorizedTable_ON = cell2table(categorizedFreezes_ON, 'VariableNames', {'FreezeType', 'StartTime', 'EndTime'});

% Display categorized freezing responses for OFF_shock
disp('Categorized Freezing Responses for OFF_shock:');
disp(categorizedTable_OFF);

% Display categorized freezing responses for ON_shock
disp('Categorized Freezing Responses for ON_shock:');
disp(categorizedTable_ON);

% Display total durations for OFF_shock
disp('Total Durations for OFF_shock:');
disp(totalDuration_OFF);

% Display total durations for ON_shock
disp('Total Durations for ON_shock:');
disp(totalDuration_ON);

%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% Define shock periods (in seconds)
OFF_shock_periods = [300, 350, 390, 420, 480, 530, 560, 600, 660, 690];
ON_shock_periods = [730, 780, 820, 880, 910, 950, 990, 1050, 1100, 1160];

N=length(Spk_data);

% Initialize counts and total observation time for OFF and ON conditions
% Initialize results structure
results_OFF = struct('NeuronID', [], 'TotalSpikes_OFF', [], 'FreezingSpikes_OFF', [], ...
                     'NonFreezingSpikes_OFF', [], 'ShockFreezingSpikes_OFF', [], ...
                     'PostShockFreezingSpikes_OFF', [], 'OutOfBoundFreezingSpikes_OFF', [], ...
                     'SpikeTimestamps_OFF', []);

results_ON = struct('NeuronID', [], 'TotalSpikes_ON', [], 'FreezingSpikes_ON', [], ...
                    'NonFreezingSpikes_ON', [], 'ShockFreezingSpikes_ON', [], ...
                    'PostShockFreezingSpikes_ON', [], 'OutOfBoundFreezingSpikes_ON', [],  ...
                    'SpikeTimestamps_ON', []);

% Loop through each neuron
for neuronIndex = 1:N
    spikeTimestamps = Spk_data(neuronIndex).Spikes; % Get spikes for the current neuron
    
    % Create off_spikeTimestamps based on the range
    off_spikeTimestamps = spikeTimestamps(spikeTimestamps >= (OFF_shock_periods(1) - 5) & ...
                                           spikeTimestamps <= (OFF_shock_periods(end) + 5));
    
    % Calculate TotalSpikes_OFF
    TotalSpikes_OFF = length(off_spikeTimestamps);

    % Initialize counts for OFF condition
    off_freezingSpikeCount = 0;
    off_nonFreezingSpikeCount = 0;
    off_shockFreezingSpikeCount = 0;
    off_postShockFreezingSpikeCount = 0;
    off_outOfBoundFreezingSpikeCount = 0;

    %TimeDuration_OFF = OFF_shock_periods(end)+5-OFF_shock_periods(1)-5; % calculate the total time duration of OFF condition



    % Initialize durations for OFF condition
    shockFreezingDuration_OFF = 0;
    postShockFreezingDuration_OFF = 0;
    totalFreezingDuration_OFF = 0;
    outOfBoundFreezingDuration_OFF = 0;


    %calculate the time duration in freeze conditions
     for i = 1:size(categorizedFreezes_OFF, 1)
        category = categorizedFreezes_OFF{i, 1};
        categoryStart = categorizedFreezes_OFF{i, 2};
        categoryEnd = categorizedFreezes_OFF{i, 3};

        % Calculate duration based on the category
        duration = categoryEnd - categoryStart;
        
        switch category
            case 'shock_freeze'
                shockFreezingDuration_OFF = shockFreezingDuration_OFF + duration;
            case 'post_shock_freeze'
                postShockFreezingDuration_OFF = postShockFreezingDuration_OFF + duration;
            case 'out_of_bound_freeze'
                outOfBoundFreezingDuration_OFF = outOfBoundFreezingDuration_OFF + duration;
        end

       % totalFreezingDuration_OFF = totalFreezingDuration_OFF + duration;
    end


    % Count spikes in OFF categorized freeze responses
    for i = 1:size(categorizedFreezes_OFF, 1)
        category = categorizedFreezes_OFF{i, 1};
        categoryStart = categorizedFreezes_OFF{i, 2};
        categoryEnd = categorizedFreezes_OFF{i, 3};

         duration = categoryEnd - categoryStart;


        % Count spikes that fall within the category time range
        for spikeIndex = 1:length(off_spikeTimestamps)
            spikeTime = off_spikeTimestamps(spikeIndex);
            if spikeTime >= categoryStart && spikeTime < categoryEnd
                switch category
                    case 'shock_freeze'
                        off_shockFreezingSpikeCount = off_shockFreezingSpikeCount + 1;
                        off_freezingSpikeCount = off_freezingSpikeCount + 1;
                    case 'post_shock_freeze'
                        off_postShockFreezingSpikeCount = off_postShockFreezingSpikeCount + 1;
                        off_freezingSpikeCount = off_freezingSpikeCount + 1;
                    case 'out_of_bound_freeze'
                        off_outOfBoundFreezingSpikeCount = off_outOfBoundFreezingSpikeCount + 1;
                         off_freezingSpikeCount = off_freezingSpikeCount + 1;
                end
            end
        end
    end

    % Count non-freezing spikes for OFF condition
    off_nonFreezingSpikeCount = TotalSpikes_OFF - off_freezingSpikeCount;

    % Store results for OFF condition

    results_OFF(neuronIndex).NeuronID = neuronIndex;
    results_OFF(neuronIndex).TotalSpikes_OFF = TotalSpikes_OFF;
    results_OFF(neuronIndex).FreezingSpikes_OFF = off_freezingSpikeCount;
    results_OFF(neuronIndex).NonFreezingSpikes_OFF = off_nonFreezingSpikeCount;
    results_OFF(neuronIndex).ShockFreezingSpikes_OFF = off_shockFreezingSpikeCount;
    results_OFF(neuronIndex).PostShockFreezingSpikes_OFF = off_postShockFreezingSpikeCount;
    results_OFF(neuronIndex).OutOfBoundFreezingSpikes_OFF = off_outOfBoundFreezingSpikeCount;
    results_OFF(neuronIndex).SpikeTimestamps_OFF = off_spikeTimestamps;
    %results_OFF(neuronIndex).TotalFreezingSpikeRate = off_freezingSpikeCount/ totalFreezingDuration_OFF;
    results_OFF(neuronIndex).ShockFreezingSpikeRate = off_shockFreezingSpikeCount/shockFreezingDuration_OFF;
    results_OFF(neuronIndex).PostShockFreezingSpikeRate = off_postShockFreezingSpikeCount/postShockFreezingDuration_OFF;
    results_OFF(neuronIndex).outOfBoundFreezingSpikeRate =  off_outOfBoundFreezingSpikeCount/outOfBoundFreezingDuration_OFF;
    results_OFF(neuronIndex).TotalSpikes= length(spikeTimestamps); 

%



    % Create on_spikeTimestamps based on the range
    on_spikeTimestamps = spikeTimestamps(spikeTimestamps >= (ON_shock_periods(1) - 5) & ...
                                           spikeTimestamps <= (ON_shock_periods(end) + 5));

    % Calculate TotalSpikes_ON
    TotalSpikes_ON = length(on_spikeTimestamps);
    
  TimeDuration_ON = ON_shock_periods(end)+5-ON_shock_periods(1)-5; % calculate the total time duration of OFF condition



 % Initialize durations for OFF condition
    shockFreezingDuration_ON = 0;
    postShockFreezingDuration_ON = 0;
   % totalFreezingDuration_ON = 0;
    outOfBoundFreezingDuration_ON = 0;


    %calculate the time duration in freeze conditions
     for i = 1:size(categorizedFreezes_ON, 1)
        category = categorizedFreezes_ON{i, 1};
        categoryStart = categorizedFreezes_ON{i, 2};
        categoryEnd = categorizedFreezes_ON{i, 3};

        % Calculate duration based on the category
        duration = categoryEnd - categoryStart;
        
        switch category
            case 'shock_freeze'
                shockFreezingDuration_ON = shockFreezingDuration_ON + duration;
            case 'post_shock_freeze'
                postShockFreezingDuration_ON = postShockFreezingDuration_ON + duration;
            case 'out_of_bound_freeze'
                outOfBoundFreezingDuration_ON = outOfBoundFreezingDuration_ON + duration;
        end

      %  totalFreezingDuration_ON = totalFreezingDuration_ON + duration;
    end


    % Initialize counts for ON condition
    on_freezingSpikeCount = 0;
    on_nonFreezingSpikeCount = 0;
    on_shockFreezingSpikeCount = 0;
    on_postShockFreezingSpikeCount = 0;
    on_outOfBoundFreezingSpikeCount = 0;

    % Count spikes in ON categorized freeze responses
    for i = 1:size(categorizedFreezes_ON, 1)
        category = categorizedFreezes_ON{i, 1};
        categoryStart = categorizedFreezes_ON{i, 2};
        categoryEnd = categorizedFreezes_ON{i, 3};

        % Count spikes that fall within the category time range
        for spikeIndex = 1:length(on_spikeTimestamps)
            spikeTime = on_spikeTimestamps(spikeIndex);
            if spikeTime >= categoryStart && spikeTime < categoryEnd
                switch category
                    case 'shock_freeze'
                        on_shockFreezingSpikeCount = on_shockFreezingSpikeCount + 1;
                        on_freezingSpikeCount = on_freezingSpikeCount + 1;
                    case 'post_shock_freeze'
                        on_postShockFreezingSpikeCount = on_postShockFreezingSpikeCount + 1;
                        on_freezingSpikeCount = on_freezingSpikeCount + 1;
                    case 'out_of_bound_freeze'
                        on_outOfBoundFreezingSpikeCount = on_outOfBoundFreezingSpikeCount + 1;
                        on_freezingSpikeCount = on_freezingSpikeCount + 1;
                end
            end
        end
    end

    % Count non-freezing spikes for ON condition
    on_nonFreezingSpikeCount = TotalSpikes_ON - on_freezingSpikeCount;

    % Store results for ON condition
    results_ON(neuronIndex).NeuronID = neuronIndex;
    results_ON(neuronIndex).TotalSpikes_ON = TotalSpikes_ON;
    results_ON(neuronIndex).FreezingSpikes_ON = on_freezingSpikeCount;
    results_ON(neuronIndex).NonFreezingSpikes_ON = on_nonFreezingSpikeCount;
    results_ON(neuronIndex).ShockFreezingSpikes_ON = on_shockFreezingSpikeCount;
    results_ON(neuronIndex).PostShockFreezingSpikes_ON = on_postShockFreezingSpikeCount;
    results_ON(neuronIndex).OutOfBoundFreezingSpikes_ON = on_outOfBoundFreezingSpikeCount;
    results_ON(neuronIndex).SpikeTimestamps_ON = on_spikeTimestamps;
   % results_ON(neuronIndex).TotalFreezingSpikeRate = on_freezingSpikeCount/ totalFreezingDuration_ON;
    results_ON(neuronIndex).ShockFreezingSpikeRate = on_shockFreezingSpikeCount/shockFreezingDuration_ON;
    results_ON(neuronIndex).PostShockFreezingSpikeRate = on_postShockFreezingSpikeCount/postShockFreezingDuration_ON;
    results_ON(neuronIndex).outOfBoundFreezingSpikeRate =  on_outOfBoundFreezingSpikeCount/outOfBoundFreezingDuration_ON;
    results_ON(neuronIndex).TotalSpikes= length(spikeTimestamps); 

end

% Display results for OFF condition
disp('OFF Condition Neuron Spike Counts and Spike Timestamps:');
disp(struct2table(results_OFF));

% Display results for ON condition
disp('ON Condition Neuron Spike Counts and Spike Timestamps:');
disp(struct2table(results_ON));


%%

% Identify the habituation period freezing
% Plot freezing responses
  habituation_freezing=0;

for i = 1:size(freezingResponses, 1)
   if freezingResponses{i,2}<300

    startTime = freezingResponses{i, 2}; % Freezing start time
    endTime = freezingResponses{i, 3};   % Freezing end time
    duration=endTime-startTime;
    habituation_freezing=habituation_freezing+duration;
   end
end







N=length(Spk_data); 

for neuronIndex = 2:N
  spikeTimestamps = Spk_data(neuronIndex).Spikes;
  habituation_freezingSpikeCount = 0;
  habituation_nonFreezingSpikeCount = 0;
 
for i = 1:size(freezingResponses, 1)   
 if freezingResponses{i,2}<300

    startTime = freezingResponses{i, 2}; % Freezing start time
    endTime = freezingResponses{i, 3};   % Freezing end time
    count=sum(  spikeTimestamps  >= startTime &   spikeTimestamps < endTime)
  habituation_freezingSpikeCount=  habituation_freezingSpikeCount+count; clear count
 end

end
 habituation_nonfreezingSpikeCount=sum(spikeTimestamps<300)-habituation_freezingSpikeCount;
 habituation(neuronIndex).freezingSpikeCount= habituation_freezingSpikeCount;
 habituation(neuronIndex).freezingSpikeRate= habituation_freezingSpikeCount/ habituation_freezing;
habituation(neuronIndex).nonfreezingSpikeCount= habituation_nonfreezingSpikeCount;
habituation(neuronIndex).nonfreezingSpikeRate= habituation_nonfreezingSpikeCount/(300- habituation_freezing);
end


%%
for neuronIndex =2:N
 spikeTimestamps = Spk_data(neuronIndex).Spikes;

%Create the response structure of Freezing in the entire timeline of shock experiment
% Create a figure for the timeline
f= fullfig;
subplot(4,1,1)
hold on

% Plot OFF_shock periods
for j = 1:length(OFF_shock)
    fill([OFF_shock(j), OFF_shock(j) + 2, OFF_shock(j) + 2, OFF_shock(j)], ...
         [0, 0, 1, 1], 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
end
hold on
% Plot ON_shock periods
for j = 1:length(ON_shock)
    fill([ON_shock(j), ON_shock(j) + 2, ON_shock(j) + 2, ON_shock(j)], ...
         [0, 0, 1, 1], 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
end

% Plot freezing responses
for i = 1:size(freezingResponses, 1)
    startTime = freezingResponses{i, 2}; % Freezing start time
    endTime = freezingResponses{i, 3};   % Freezing end time
    plot([startTime, endTime], [0.5, 0.5], 'k', 'LineWidth', 100); % Black line for freezing
end

 trial{1}=0.66*ones(size(spikeTimestamps))
plot_raster(spikeTimestamps,cell2mat(trial(:)),0.2,'magenta')

% Set the axis limits and labels
xlim([0, 1200]);
ylim([-0.1, 1.1]);
xlabel('Time (seconds)');
ylabel('Events');
title('Timeline of OFFshock, ONshock, and Freezing Responses');
%legend({'OFFshock', 'ONshock', 'Freezing Responses'}, 'Location', 'NorthEast');
grid on;
grid minor
%hold off;
%set(gca, 'units', 'normalized','position',[0.13 0.11 0.77 0.81] )

%bar([results_OFF(neuronIndex).TotalFreezingSpikeRate results_ON(neuronIndex).TotalFreezingSpikeRate] );



  % Observed data
 n1 = results_OFF(neuronIndex).FreezingSpikes_OFF; N1 = results_OFF(neuronIndex).NonFreezingSpikes_OFF+n1;
 n2 = results_ON(neuronIndex).FreezingSpikes_ON; N2 =results_ON(neuronIndex).NonFreezingSpikes_ON+n2;
 x1 = [repmat('a',N1,1); repmat('b',N2,1)];
 x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
% x2 = [repmat(1,n1,1); repmat(2,N1,1); repmat(1,n2,1); repmat(2,N2,1)];
 if length(x1)-length(x2)==0
 [tbl,chi2stat,pval] = crosstab(x1,x2)
cnames = {'Freezing Spike Count','Non Freezing Spike Count'};
rnames = {'Observed (OFF)','Expected (ON)'};

% figure
% h = subplot(2,2,4); 

% h.Visible='off';
%subplot(2,1,1),plot(3)
% pos = get(subplot(2,1,2),'position');
% delete(subplot(2,1,2))
% %set(t,'units','normalized')
% set(t,'position',pos)

%f=figure(4)
t = uitable('Data',tbl,...
            'ColumnName',cnames,... 
            'RowName',rnames,...
            'ColumnWidth',{150})
% subplot(2,1,1),plot(3)
% 
% pos = get(subplot(2,1,2),'position');
% title('Chi Square table')
% %xlabel('xlabel')
% %ylabel('ylabel')
% set(subplot(2,1,2),'yTick',[])
% set(subplot(2,1,2),'xTick',[])
% set(t,'units','normalized')
% set(t,'position',pos)
 set(t,'ColumnName',{'Freezing Spike Count','Non Freezing Spike Count'})
set(t,'RowName',{'Observed (OFF)','Expected (ON)'})
t.Position = [30 950 t.Extent(3) t.Extent(4)];
set(gcf,'units','normalized','outerposition',[0 0 1 1])
txt= strcat('ChiSquareStat= ',num2str(chi2stat),' pvalue= ',num2str(pval)) 
annotation('textbox', [0.15, 0.74,0, 0], 'string', txt,'Color','magenta','FontSize',20)

if pval<0.05
    txt1='significant'
annotation('textbox', [0.25, 0.7,0, 0], 'string', txt1,'Color','red','FontSize',20)
else
    txt1='nonsignificant'
annotation('textbox', [0.25, 0.7,0, 0], 'string', txt1,'Color','blue','FontSize',20)
end
subplot(4,1,2)
off_freeze=[results_OFF(neuronIndex).FreezingSpikes_OFF results_OFF(neuronIndex).NonFreezingSpikes_OFF];
on_freeze= [results_ON(neuronIndex).FreezingSpikes_ON results_ON(neuronIndex).NonFreezingSpikes_ON];

FreezeComp=[off_freeze; on_freeze];
p2=bar(FreezeComp,'stacked');
title('Stacked Plot of Freezing and Non Freezing Spike Count in OFF and ON Condition')
legend('Freeze Spike Counts', 'Non Freeze Spike Counts')
xticklabels({'OFF Condition','ON Condition'})
subplot(4,1,3)
shock_freeze=[results_OFF(neuronIndex).ShockFreezingSpikes_OFF results_ON(neuronIndex).ShockFreezingSpikes_ON];
p2=bar(shock_freeze);
title('Bar graph of Shock Freeze Counts in OFF and ON Condition')
xticklabels({'OFF Condition','ON Condition'})
subplot(4,1,4)
postshock_freeze=[results_OFF(neuronIndex).PostShockFreezingSpikes_OFF results_ON(neuronIndex).PostShockFreezingSpikes_ON];

p2=bar(postshock_freeze);
title('Bar graph of Post Shock Freeze Counts in OFF and ON Condition')

xticklabels({'OFF Condition','ON Condition'})
sgtitle(strcat('NeuronID#',num2str(neuronIndex)) )

 end
end
%%
figure
hold on
for k=1:length(freezingResponses)
startX=freezingResponses{k,2}
endX=freezingResponses{k,3}
x = [startX endX endX startX];
y = [N N 2 2];
fill(x,y,'cyan')
end
hold on
     for neuronIndex=2:N

   spikeTimestamps = Spk_data(neuronIndex).Spikes;
 % plot_raster(spikeTimestamps,xlim([0 1]))

 trial{1}=neuronIndex*ones(size(spikeTimestamps));
plot_raster(spikeTimestamps,cell2mat(trial(:)),1,'k')
xline(OFF_shock,'b')
xline(ON_shock,'r')

     end
