%%%CODE SAMPLE. SCRIPT VISUALIZES SIZE OF TOPOLOGICAL OBJECTS IN RAT EEG AT DIFFERENT STATES OF CONSCIOUSNESS



bettispath = 'D:\Joe\TDA_for_rat_data_summer_2024\PyCliqueTop_2023-main\TDA_results\';
figurepath = 'D:\Joe\TDA_for_rat_data_summer_2024\PyCliqueTop_2023-main\areas_under_curve_9rat\';
addpath(bettispath);
files = ls(fullfile(bettispath, '*.mat'));
files = cellstr(files);
bandname2={'DELTA', 'THETA', 'ALPHA', 'BETA', 'GAMMA1', 'GAMMA2','GAMMA3', 'BROAD'};
band_freq = {'1-4hz','4-10hz','10-15hz','15-25hz','25-55hz','85-125hz','125-155hz','1-155hz'};
states = {'N_W', 'S_1', 'S_2', 'RoC'};
verbessStates = {'Normal Wake', 'Recovery', 'Sevo Epoch 1', 'Sevo Epoch 2'};
state_colors = {'Red', 'Green', 'Blue', 'magenta'};
ylims={30,25,20,18,12,4, 1234523613613};



for bp = 1:8   
bp_searchtext = sprintf('%.2f', bp);
    betFiles = files(contains(files,bp_searchtext));
    betFiles = reshape(betFiles, numel(betFiles)/numel(states), []);

    if bp ==1 %%find min edge density for interp1

        findMaxED = [];

        for edge = 1:size(betFiles,1)


            load(betFiles{edge,1})
            disp(betFiles{edge,1})
            avgBettiCurves = mean(BettiCurves,1);
            sqBettiCurves = squeeze(avgBettiCurves);
            EE = Betti2EulerEntropy(sqBettiCurves);
            findMaxED = [findMaxED; numel(EE)];

        end
        maxEdgeDensity = linspace(0,1,max(findMaxED));
    end 
    
    %found max edge density resolution that everything needs to be
    %interpolated to....


    bettisAcrossSets= zeros(numel(maxEdgeDensity), size(sqBettiCurves,2),size(betFiles,1)); %%edge density x bnum x dataset x num(band positions)
    for m = 1:size(betFiles,2) %%for each state
        disp(['state' num2str(m)])
        for p = 1:size(betFiles,1) %%for each data set at each state (length 9 for 9 datasets)
            load(betFiles{p,m})
            disp(betFiles{p,m})
            avgBettiCurves = mean(BettiCurves,1);
            sqBettiCurves = squeeze(avgBettiCurves);
            EE = Betti2EulerEntropy(sqBettiCurves);
            interpBettiCurves = interp1(edgeDensity, sqBettiCurves, maxEdgeDensity);
            interpEE(:,p) = interp1(edgeDensity, EE, maxEdgeDensity);

            for bet =1:size(sqBettiCurves,2)
                bettisAcrossSets(:,bet,p) = interpBettiCurves(:,bet);
            end

        end   %%all sets for one state are loaded and betti numbers are averaged
        % avgEEAcrossSets(:,m) = mean(interpEE,2);
        avgBettisAcrossSets(:,:,m,bp) = mean(bettisAcrossSets,3); %%avg figure for bettis 0-5 of all sets, third dimension is across states for the figure
    end
end

%avgBettisAcrossSets is 4-d matrix which is (edge density, avg betti Number across sets,
%conscious state, band position. We use this to make 3d plots
%for a plot,  (second dimension of
%avgBettisAcrossSets = betti number), (third dimension of avgBettisAcrossSets = conscious state) (4th dimension of
%avgBettisAcrossSets = band position), , and
%the 3rd dimension of this plot will be the areas under the curves. 


for betti = 1:size(avgBettisAcrossSets,2) %second dimension of avgBettisAcrossSets
    figure(); 
for band_position = 1:size(avgBettisAcrossSets,4) %fourth dimension of avgBettisAcrossSets
for state = 1:size(avgBettisAcrossSets,3) %third dimension of avgBettisAcrossSets
auc_for_conscious_states(state, band_position) = sum(avgBettisAcrossSets(:,betti,state,band_position));  %%going to be anxious to make sure you know what state you're working with 
end
end
bar3(auc_for_conscious_states);title(['Areas Under Curves for BettiNumber' num2str(betti-1)]);
yticklabels(verbessStates);
xticklabels(band_freq);
colormap jet
set(gca,'FontName', 'Times New Roman')
saveas(gca, (sprintf('%sAreas_under_curve_betti=%.2f.png', figurepath, betti)))
end
