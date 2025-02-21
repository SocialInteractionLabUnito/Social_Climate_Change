% ------------------------------------------------------------------------%
% Here we take pre-processed EDA data, behavioral data, subjective
% measures, stories
% We remove outlier cases from EDA data, behavioral data, subjective 
% measures (Valence and Future). 
% We test normality assumptions for EDA data, behavioral data, subjective 
% measures (Valence and Future).
% We plot data
% ------------------------------------------------------------------------%

clc
clear

% Path for data
rootdir = 'C:\Users\reven\OneDrive\Desktop\Giulia Romano Cappi\Go green\Github';

%% EDA

% Load data
eda = readtable(fullfile(rootdir,'EDA.xlsx'));
% Keep Phasic Data to plot
eda = eda(:,["Subj","Bella_Zscore_phasic","Brutta_Zscore_phasic"]);

% Remove outlier cases (within Clean and within Polluted)
isup_bella = eda.Bella_Zscore_phasic > mean(eda.Bella_Zscore_phasic) + 2*(std(eda.Bella_Zscore_phasic));
isdown_bella = eda.Bella_Zscore_phasic < mean(eda.Bella_Zscore_phasic) - 2*(std(eda.Bella_Zscore_phasic));
isup_brutta = eda.Brutta_Zscore_phasic > mean(eda.Brutta_Zscore_phasic) + 2*(std(eda.Brutta_Zscore_phasic));
isdown_brutta = eda.Brutta_Zscore_phasic < mean(eda.Brutta_Zscore_phasic) - 2*(std(eda.Brutta_Zscore_phasic));

eda(isup_bella | isdown_bella | isup_brutta | isdown_brutta,:) = [];

% Normality test for Clean
[H, pValue, SWstatistic] = swtest(eda.Bella_Zscore_phasic, 0.05)
% Normality test for Clean
[H, pValue, SWstatistic] = swtest(eda.Brutta_Zscore_phasic, 0.05)
% T-test
[h,p,ci,stats] = ttest(eda.Bella_Zscore_phasic,eda.Brutta_Zscore_phasic)

% plot
% violin
figure
ydata = [eda.Bella_Zscore_phasic,eda.Brutta_Zscore_phasic];
xdata = categorical(["Bella","Brutta"]);
violinplot(xdata,ydata)
hold on
yline(median(eda.Bella_Zscore_phasic),'b')
yline(mean(eda.Bella_Zscore_phasic),'g')
yline(median(eda.Brutta_Zscore_phasic),'r')
yline(mean(eda.Brutta_Zscore_phasic),'c')
ylabel("EDA Mean [Zscore]")
title("EDA")
swarmchart(xdata,[eda.Bella_Zscore_phasic,eda.Brutta_Zscore_phasic],"black","filled")
for i = 1:size(eda,1)
    currdata = horzcat(table2array(eda(i,"Bella_Zscore_phasic")),table2array(eda(i,"Brutta_Zscore_phasic")));
    plot(currdata,"black")
end
legend("","","Median Bella","Mean Bella","Median Brutta","Mean Brutta","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","")


%% N photos

% Load data
photos = readtable(fullfile(rootdir,'AllBehavioral.xlsx'));
% Keep only pictures
photos = photos(:,["Subj","N_Photo_good","N_Photo_bad"]);

% Remove outlier cases (within Clean and within Polluted)
isup_bella = photos.N_Photo_good > mean(photos.N_Photo_good) + 2*(std(photos.N_Photo_good));
isdown_bella = photos.N_Photo_good < mean(photos.N_Photo_good) - 2*(std(photos.N_Photo_good));
isup_brutta = photos.N_Photo_bad > mean(photos.N_Photo_bad) + 2*(std(photos.N_Photo_bad));
isdown_brutta = photos.N_Photo_bad < mean(photos.N_Photo_bad) - 2*(std(photos.N_Photo_bad));

photos(isup_bella | isdown_bella | isup_brutta | isdown_brutta,:) = [];

% Normality test for Clean
[H, pValue, SWstatistic] = swtest(photos.N_Photo_good, 0.05)
% Normality test for Clean
[H, pValue, SWstatistic] = swtest(photos.N_Photo_bad, 0.05)
% Wilcoxon signed rank test
[p,h,stats] = signrank(photos.N_Photo_good,photos.N_Photo_bad)

% plot
% bar
data = [mean(photos.N_Photo_good),mean(photos.N_Photo_bad)];
err = [(std(photos.N_Photo_good,'omitnan') / sqrt(length(photos.N_Photo_good))),(std(photos.N_Photo_bad,'omitnan') / sqrt(length(photos.N_Photo_bad)))];

figure
x = categorical({'Bella','Brutta'});
x = reordercats(x,{'Bella','Brutta'});
bar(x,data)
hold on
er = errorbar(x,data,err);
er.Color = [0 0 0];
er.LineStyle = 'none';
ylabel("N Photos")
title("Photos")
ylim([3,6])


%% N movements

% Load data
mov = readtable(fullfile(rootdir,'AllBehavioral.xlsx'));
% Keep only movements
mov = mov(:,["Subj","N_Teleportation_Good","N_Teleportation_Bad"]);

% Remove outlier cases (within Clean and within Polluted)
isup_bella = mov.N_Teleportation_Good > mean(mov.N_Teleportation_Good) + 2*(std(mov.N_Teleportation_Good));
isdown_bella = mov.N_Teleportation_Good < mean(mov.N_Teleportation_Good) - 2*(std(mov.N_Teleportation_Good));
isup_brutta = mov.N_Teleportation_Bad > mean(mov.N_Teleportation_Bad) + 2*(std(mov.N_Teleportation_Bad));
isdown_brutta = mov.N_Teleportation_Bad < mean(mov.N_Teleportation_Bad) - 2*(std(mov.N_Teleportation_Bad));

mov(isup_bella | isdown_bella | isup_brutta | isdown_brutta,:) = [];

% Normality test for Clean
[H, pValue, SWstatistic] = swtest(mov.N_Teleportation_Good, 0.05)
% Normality test for Clean
[H, pValue, SWstatistic] = swtest(mov.N_Teleportation_Bad, 0.05)
% Wilcoxon signed rank test
[p,h,stats] = signrank(mov.N_Teleportation_Good,mov.N_Teleportation_Bad)

% plot
% bar
data = [mean(mov.N_Teleportation_Good),mean(mov.N_Teleportation_Bad)];
err = [(std(mov.N_Teleportation_Good,'omitnan') / sqrt(length(mov.N_Teleportation_Good))),(std(mov.N_Teleportation_Bad,'omitnan') / sqrt(length(mov.N_Teleportation_Bad)))];

figure
x = categorical({'Bella','Brutta'});
x = reordercats(x,{'Bella','Brutta'});
bar(x,data)
hold on
er = errorbar(x,data,err);
er.Color = [0 0 0];
er.LineStyle = 'none';
ylabel("N Movements")
title("Movements")
ylim([40,70])

% Let's control for possible interaction effects between EDA and Movements
% Set data for model

% EDA
neweda = table();
neweda.Subj = repmat(eda.Subj,2,1);
neweda.Environment = vertcat(zeros(1,size(eda,1))',ones(1,size(eda,1))');
neweda.EDA = vertcat(eda.Bella_Zscore_phasic,eda.Brutta_Zscore_phasic);

% Movements
newmov = table();
newmov.Subj = repmat(mov.Subj,2,1);
newmov.Environment = vertcat(zeros(1,size(mov,1))',ones(1,size(mov,1))');
newmov.Movements = vertcat(mov.N_Teleportation_Good,mov.N_Teleportation_Bad);

% Let's create a table with the same subjects for EDA and for Movements
sogg_eda = unique(neweda.Subj);
sogg_tp = unique(newmov.Subj);
sogg_shared = intersect(sogg_eda,sogg_tp);

if ~isempty(sogg_shared)
    neweda = neweda(ismember(neweda.Subj, sogg_shared), :);
    newmov = newmov(ismember(newmov.Subj, sogg_shared), :);
end

newdata = horzcat(neweda(:,:),newmov(:,"Movements"));

% Let's first run a model with: 
% EDA as dependent measure,
% Environment (Clean/Polluted) and N movements as fixed effects,
% Subject as random effect
 
mdl = fitlme(newdata, 'EDA ~ Environment*Movements + (1|Subj)')

% split data for condition
eda_beh_bella = newdata(newdata.Environment == 0,:);
eda_beh_brutta = newdata(newdata.Environment == 1,:);

% Re-test for Clean only
mdl = fitlme(eda_beh_bella, 'EDA ~ Movements + (1|Subj)') 
% Re-test for Polluted only
mdl = fitlme(eda_beh_brutta, 'EDA ~ Movements + (1|Subj)') 


%% Valence

% Load data
valence = readtable(fullfile(rootdir,'SubjectiveMeasures.xlsx'));
% Keep only valence
valence = valence(:,["Subj","Emo_NonPolluted","Emo_Polluted"]);

% Remove nans if present
valence = valence(~any(isnan(valence{:, 2:3}), 2), :);

% Remove outlier cases (within Clean and within Polluted)
isup_bella = valence.Emo_NonPolluted > mean(valence.Emo_NonPolluted) + 2*(std(valence.Emo_NonPolluted));
isdown_bella = valence.Emo_NonPolluted < mean(valence.Emo_NonPolluted) - 2*(std(valence.Emo_NonPolluted));
isup_brutta = valence.Emo_Polluted > mean(valence.Emo_Polluted) + 2*(std(valence.Emo_Polluted));
isdown_brutta = valence.Emo_Polluted < mean(valence.Emo_Polluted) - 2*(std(valence.Emo_Polluted));

valence(isup_bella | isdown_bella | isup_brutta | isdown_brutta,:) = [];

% Normality test for Clean
[H, pValue, SWstatistic] = swtest(mov.N_Teleportation_Good, 0.05)
% Normality test for Clean
[H, pValue, SWstatistic] = swtest(mov.N_Teleportation_Bad, 0.05)
% Wilcoxon signed rank test
[p,h,stats] = signrank(valence.Emo_NonPolluted,valence.Emo_Polluted)


%% Future

% Load data
future = readtable(fullfile(rootdir,'SubjectiveMeasures.xlsx'));
% Keep only future
future = future(:,["Subj","Futuro_Reale_NonPolluted","Futuro_Reale_Polluted"]);

% Remove nans if present
future = future(~any(isnan(future{:, 2:3}), 2), :);

% Remove outlier cases (within Clean and within Polluted)
isup_bella = future.Futuro_Reale_NonPolluted > mean(future.Futuro_Reale_NonPolluted) + 2*(std(future.Futuro_Reale_NonPolluted));
isdown_bella = future.Futuro_Reale_NonPolluted < mean(future.Futuro_Reale_NonPolluted) - 2*(std(future.Futuro_Reale_NonPolluted));
isup_brutta = future.Futuro_Reale_Polluted > mean(future.Futuro_Reale_Polluted) + 2*(std(future.Futuro_Reale_Polluted));
isdown_brutta = future.Futuro_Reale_Polluted < mean(future.Futuro_Reale_Polluted) - 2*(std(future.Futuro_Reale_Polluted));

future(isup_bella | isdown_bella | isup_brutta | isdown_brutta,:) = [];

% Normality test for Clean
[H, pValue, SWstatistic] = swtest(future.Futuro_Reale_NonPolluted, 0.05)
% Normality test for Clean
[H, pValue, SWstatistic] = swtest(future.Futuro_Reale_Polluted, 0.05)
% Wilcoxon signed rank test
[p,h,stats] = signrank(future.Futuro_Reale_NonPolluted,future.Futuro_Reale_Polluted)

% Plot Valence and Future together

% Valence
data_v = [mean(valence.Emo_NonPolluted),mean(valence.Emo_Polluted)];
err_v = [(std(valence.Emo_NonPolluted,'omitnan') / sqrt(length(valence.Emo_NonPolluted))),(std(valence.Emo_Polluted,'omitnan') / sqrt(length(valence.Emo_Polluted)))];
% Future
data_f = [mean(future.Futuro_Reale_NonPolluted),mean(future.Futuro_Reale_Polluted)];
err_f = [(std(future.Futuro_Reale_NonPolluted,'omitnan') / sqrt(length(future.Futuro_Reale_NonPolluted))),(std(future.Futuro_Reale_Polluted,'omitnan') / sqrt(length(future.Futuro_Reale_Polluted)))];

figure
x = categorical({'Bella','Brutta'});
x = reordercats(x,{'Bella','Brutta'});
yyaxis left
plot(x,data_v)
hold on
er = errorbar(x,data_v,err_v);
er.Color = [0 0 0];
er.LineStyle = 'none';
ylabel("Valence")
ylim([1,7])
yyaxis right
plot(x,data_f)
er = errorbar(x,data_f,err_v);
er.Color = [0 0 0];
er.LineStyle = 'none';
ylabel("Future")
ylim([2,6])
title("Subjective measures")

%% Subjective Ratings

% Vision
% Load data
data = readtable(fullfile(rootdir,'SubjectiveMeasures.xlsx'));
% Keep only vision
vision = data(:,["Subj","Qualita_Visione"]);
% Remove nans if present
vision(isnan(vision.Qualita_Visione),:) = [];
% Extract mean and standard deviation
mean_vision = mean(vision.Qualita_Visione);
sd_vision = std(vision.Qualita_Visione);

% Sound
% keep only sounds
sounds = data(:,["Subj","Qualita_Suoni"]);
% Remove nans if present
sounds(isnan(sounds.Qualita_Suoni),:) = [];
% Extract mean and standard deviation
mean_sound = mean(sounds.Qualita_Suoni);
sd_sound = std(sounds.Qualita_Suoni);

% Movement
% keep only movement
movement = data(:,["Subj","Qualita_Movimento"]);
% Remove nans if present
movement(isnan(movement.Qualita_Movimento),:) = [];
% Extract mean and standard deviation
mean_mov = mean(movement.Qualita_Movimento);
sd_mov = std(movement.Qualita_Movimento);

% Photos
% keep only photos
photoq = data(:,["Subj","Qualita_Foto"]);
% Remove nans if present
photoq(isnan(photoq.Qualita_Foto),:) = [];
% Extract mean and standard deviation
mean_phq = mean(photoq.Qualita_Foto);
sd_phq = std(photoq.Qualita_Foto);

% Realism
% keep only realism
realism = data(:,["Subj","Realismo_Piazza"]);
% Remove nans if present
realism(isnan(realism.Realismo_Piazza),:) = [];
% Extract mean and standard deviation
mean_real = mean(realism.Realismo_Piazza);
sd_real = std(realism.Realismo_Piazza);

% Cybersyckness
% keep only photos
cysy = data(:,["Subj","Fastidio"]);
% Remove nans if present
cysy(isnan(cysy.Fastidio),:) = [];
% Extract mean and standard deviation
mean_cysy = mean(cysy.Fastidio);
sd_cysy = std(cysy.Fastidio);

% Enjoyment
% keep only photos
bore = data(:,["Subj","Noia"]);
% Remove nans if present
bore(isnan(bore.Noia),:) = [];
% Extract mean and standard deviation
mean_bore = mean(bore.Noia);
sd_bore = std(bore.Noia);

% Environmental Awareness
% keep only photos
env = data(:,["Subj","CC_ProblemaOggi"]);
% Remove nans if present
env(isnan(env.CC_ProblemaOggi),:) = [];
% Extract mean and standard deviation
mean_env = mean(env.CC_ProblemaOggi);
sd_env = std(env.CC_ProblemaOggi);

%% Stories Validation

% Load data
data = readtable(fullfile(rootdir,'Stories.xlsx'));

% Keep only anxious ratings
data_anx = table2array(data(:,contains(data.Properties.VariableNames,'anx')));

% Keep only hopeful ratings
data_hop = table2array(data(:,contains(data.Properties.VariableNames,'hop')));

% Keep only valence ratings
data_val = table2array(data(:,contains(data.Properties.VariableNames,'emo')));

% Keep only arousal ratings
data_aro = table2array(data(:,contains(data.Properties.VariableNames,'fisio')));


data_new = [mean(data_anx(:,1)), mean(data_hop(:,1)), mean(data_val(:,1)), mean(data_aro(:,1)); 
            mean(data_anx(:,2)), mean(data_hop(:,2)), mean(data_val(:,2)), mean(data_aro(:,2));
            mean(data_anx(:,3)), mean(data_hop(:,3)), mean(data_val(:,3)), mean(data_aro(:,3));
            mean(data_anx(:,4)), mean(data_hop(:,4)), mean(data_val(:,4)), mean(data_aro(:,4)); 
            mean(data_anx(:,5)), mean(data_hop(:,5)), mean(data_val(:,5)), mean(data_aro(:,5));
            mean(data_anx(:,6)), mean(data_hop(:,6)), mean(data_val(:,6)), mean(data_aro(:,6));
            mean(data_anx(:,7)), mean(data_hop(:,7)), mean(data_val(:,7)), mean(data_aro(:,7));
            mean(data_anx(:,8)), mean(data_hop(:,8)), mean(data_val(:,8)), mean(data_aro(:,8))];

stories = {'Story 1', 'Story 2', 'Story 3', 'Story 4', 'Story 5', 'Story 6', 'Story 7', 'Story 8'};
dimensions = {'Anxiety', 'Hopefulness', 'Valence', 'Arousal'};

% Plot the heatmap
figure
h = heatmap(dimensions, stories, data_new, 'ColorbarVisible', 'on');
h.Colormap = turbo;
title('Mean Ratings for Anxiety, Hopefulness, Valence, and Arousal');
xlabel('Dimension');
ylabel('Stories');



