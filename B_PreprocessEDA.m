% ------------------------------------------------------------------------%
% Here we take segmented data previously obtained from raw data (Raw data
% are available upon request). We pre-process data, zscore it, and extract
% mean values for subsequent analyses --> See C_Stats_Plots.m for stats and
% plots
% ------------------------------------------------------------------------%

clc
clear

% Set stuff
%------------------------------------------------------------------
% Load data
rootdir = 'C:\Users\reven\OneDrive\Desktop\Giulia Romano Cappi\Go green\Github';
load(fullfile(rootdir, 'AllSignals_Raw.mat'))

% Sampling rate
SR = 500;

% Options
removesoggOption = 1; % remove bad subjects
removefinalOption = 1; % remove final 90 seconds
cutBLOption = 1; % cut all baselines to 60 seconds
cutdataOption = 1; % cut data to the shortest
filterOption = 1; % filter data
ZscoreBLOption_order = 1; % zscore the signal based on the baseline

% ------------------------------------------------------------------------%

% Remove bad subjects
if removesoggOption 
    bad_subj = ["Sogg_01"];
    if ~isempty(bad_subj)
        allSignals(ismember(allSignals.Subj, bad_subj), :) = [];
    end
else
end

% Remove final 90 seconds
if removefinalOption
    % Bad condition
    bad_cond = ["BaselineFinale"];
    % Clear bad condition
    if ~isempty(bad_cond)
        allSignals(ismember(allSignals.Condition, bad_cond), :) = [];
    end
end

sogg = unique(allSignals.Subj);

% Cut all baselines to 60 seconds
if cutBLOption
    newdata = table();
    for i = [sogg]'
        currsogg = allSignals(allSignals.Subj == i,:);
        isBLBella = currsogg.Condition == "BaselineBella";
        isBLBrutta = currsogg.Condition == "BaselineBrutta";

        if size(currsogg.Raw{isBLBella},1) > (60*SR) | size(currsogg.Raw{isBLBrutta},1) > (60*SR)

            currsogg.Raw{isBLBella} = currsogg.Raw{isBLBella}(end-SR*60:end);
            currsogg.Raw{isBLBrutta} = currsogg.Raw{isBLBrutta}(end-SR*60:end);
        else
        end

        newdata = vertcat(newdata,currsogg);
    end
    allSignals = newdata;
end

% Cut data to the shortest
if cutdataOption

    sizes = table();
    for i = 1:height(allSignals)
        currline = allSignals(i,:);
        currline.size = size(currline.Raw{1},1);
        sizes = vertcat(sizes,currline);
    end

    isBellaBrutta = sizes.Condition == "Bella" | sizes.Condition == "Brutta";
    minsize = min(sizes.size(isBellaBrutta));
    NewData = table();
    for i = 1:height(allSignals)
        currline = allSignals(i,:);
        if size(currline.Raw{1,1},1) < minsize
            currline.Raw{1,1} = currline.Raw{1,1};
        elseif size(currline.Raw{1,1},1) ~= minsize
            currline.Raw{1,1} = currline.Raw{1,1}(end-minsize:end-1);
        else
        end
        NewData = vertcat(NewData,currline);
    end

    allSignals = NewData;

end

% Filter data
if filterOption 
    for i = 1:size(allSignals,1)
        currdata = allSignals.Raw{i,1};
        tx = (0:size(currdata,1)-1)'*1/SR;
        currdata = smoothdata(currdata,'movmean',0.5,'SamplePoints',tx);
        allSignals.SCL{i,1} = currdata;
        currdata = highpass(currdata,0.05,SR);
        currdata = lowpass(currdata,20,SR);
        allSignals.Phasic{i,1} = currdata;
    end
end

% Zscore the signal based on the baseline
if ZscoreBLOption_order
    allSignals_order = table();
    allSignals_nuevo = allSignals;
    for j = [sogg]'

        % Search baselines
        for m = 1:size(allSignals_nuevo,1)
            if allSignals_nuevo.Condition(m) == "BaselineBella" & allSignals_nuevo.sequenza(m) == 1
                allSignals_nuevo.Condition(m) = "BaselinePrima";
            elseif allSignals_nuevo.Condition(m) == "BaselineBrutta" & allSignals_nuevo.sequenza(m) == 1
                allSignals_nuevo.Condition(m) = "BaselineSeconda";
            elseif allSignals_nuevo.Condition(m) == "BaselineBella" & allSignals_nuevo.sequenza(m) == 2
                allSignals_nuevo.Condition(m) = "BaselineSeconda";
            elseif allSignals_nuevo.Condition(m) == "BaselineBrutta" & allSignals_nuevo.sequenza(m) == 2
                allSignals_nuevo.Condition(m) = "BaselinePrima";
            else
            end
        end

        currsogg = allSignals_nuevo(allSignals_nuevo.Subj == j,:);

        % Get indices
        isBLPrima = currsogg.Condition == "BaselinePrima";
        isBLSeconda = currsogg.Condition == "BaselineSeconda";
        isBella = currsogg.Condition == "Bella";
        isBrutta = currsogg.Condition == "Brutta";

        % Zscore phasic
        currsogg.Zscore_phasic{isBLPrima,1} = currsogg.Phasic{isBLPrima,1};
        currsogg.Zscore_phasic{isBLSeconda,1} = currsogg.Phasic{isBLSeconda,1};
        meanBLPrima = mean(currsogg.Phasic{isBLPrima,1});
        sdBLPrima = std(currsogg.Phasic{isBLPrima,1});

        currsogg.Zscore_phasic{isBella,1} = (currsogg.Phasic{isBella,1}-meanBLPrima)/sdBLPrima;
        currsogg.Zscore_phasic{isBrutta,1} = (currsogg.Phasic{isBrutta,1}-meanBLPrima)/sdBLPrima;

        % Zscore tonic
        currsogg.Zscore_tonic{isBLPrima,1} = currsogg.SCL{isBLPrima,1};
        currsogg.Zscore_tonic{isBLSeconda,1} = currsogg.SCL{isBLSeconda,1};
        meanBLPrima = mean(currsogg.SCL{isBLPrima,1});
        sdBLPrima = std(currsogg.SCL{isBLPrima,1});

        currsogg.Zscore_tonic{isBella,1} = (currsogg.SCL{isBella,1}-meanBLPrima)/sdBLPrima;
        currsogg.Zscore_tonic{isBrutta,1} = (currsogg.SCL{isBrutta,1}-meanBLPrima)/sdBLPrima;

        % Assign
        allSignals_order = vertcat(allSignals_order,currsogg);
    end
end

allSignals = horzcat(allSignals,allSignals_order.Zscore_phasic,allSignals_order.Zscore_tonic);
allSignals = renamevars(allSignals,"Var7","Zscore_phasic");
allSignals = renamevars(allSignals,"Var8","Zscore_tonic");


% Create excel file with all measures
measures_i_want = ["Raw", "SCL", "Phasic", "Zscore_phasic","Zscore_tonic"];
Results = table(); % qua ci mettiamo le misure medie (Raw, SCL, Phasic, Zscore, ZscoreSCL) in un foglio excel
for j = [sogg]'
    results = table();
    currsogg = allSignals(allSignals.Subj == j,:);
    numsogg = str2double(extractAfter(currsogg.Subj(1),'_'));
    results.Subj = currsogg.Subj(1);
    results.Sequenza = currsogg.sequenza(1);
    isBLBella = currsogg.Condition == "BaselineBella";
    isBLBrutta = currsogg.Condition == "BaselineBrutta";
    isBella = currsogg.Condition == "Bella";
    isBrutta = currsogg.Condition == "Brutta";
    fns = string(fieldnames(currsogg));
    for m = 1:numel(measures_i_want)
        myM = measures_i_want(m);
        idx = find(fns==myM);
        currVar = fns(idx);
        newVar = ['BLBella_' char(currVar)];
        results.(newVar) = mean(currsogg.(currVar){isBLBella,1}, 'omitnan');
        newVar = ['BLBrutta_' char(currVar)];
        results.(newVar) = mean(currsogg.(currVar){isBLBrutta,1}, 'omitnan');
        newVar = ['Bella_' char(currVar)];
        results.(newVar) = mean(currsogg.(currVar){isBella,1}, 'omitnan');
        newVar = ['Brutta_' char(currVar)];
        results.(newVar) = mean(currsogg.(currVar){isBrutta,1}, 'omitnan');
    end
    
    Results = vertcat(Results,results);
end

writetable(Results, [rootdir, '\EDA.xlsx']);

