%%%CODE SAMPLE. TEST FOR SIGNIFICANT CHANGES IN THE TOPOLOGICAL PATTERNS OF RAT EEG AFTER ANESTHESIA AND RECOVERY. 
%%%THIS ANALYSIS GAVE A STATISTICALLY SIGNIFICANT (AND PUBLISHABLE) RESULT




bettispath = 'D:\Joe\TDA_for_rat_data_summer_2024\PyCliqueTop_2023-main\TDA_results\';
figurepath = 'D:\Joe\TDA_for_rat_data_summer_2024\PyCliqueTop_2023-main\prelim_ttest_areas_9rat\';
addpath(bettispath);
files = ls(fullfile(bettispath, '*.mat'));
files = cellstr(files);
bandname2={'DELTA', 'THETA', 'ALPHA', 'BETA', 'GAMMA1', 'GAMMA2','GAMMA3', 'BROAD'};
band_freq = {'1-4hz','4-10hz','10-15hz','15-25hz','25-55hz','85-125hz','125-155hz','1-155hz'};
states = {'N_W', 'S_1', 'S_2', 'RoC'};
verbessStates = {'Normal Wake', 'Recovery', 'Sevo Epoch 1', 'Sevo Epoch 2'};
state_colors = {'Red', 'Green', 'Blue', 'magenta'};
ylims={30,25,20,18,12,4, 1234523613613};



for bp = 1:8 %%for only the theta band
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
    
    %found max edge density resolution to which everything must be
    %interpolated 

  
    clc;

%edge density x % bnum x dataset x number of states
bettisAcrossSets= zeros(numel(maxEdgeDensity), size(sqBettiCurves,2),size(betFiles,1),numel(states));
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
                bettisAcrossSets(:,bet,p,m) = interpBettiCurves(:,bet);
            end
        end 
    end



%%for a betti number, for a frequency band. Remember that bettisAcrossSets format is
% (edge density x % bnum x dataset x number of states)

for betti = 1:size(bettisAcrossSets,2)
figure();
ttestmats = squeeze(sum(squeeze(bettisAcrossSets(:,betti,:,:)))); 
ttestresult = zeros(size(ttestmats,2),size(ttestmats,2));
for i = 1:size(ttestmats,2)
    for j = 1:size(ttestmats,2)
[R,p]=ttest2(ttestmats(:,i),ttestmats(:,j));
if p < .05
ttestresult(i,j) = 1;
end
if p < .01
    ttestresult(i,j) = .5;
end
    end
imagesc(ttestresult); title(['ttest result for Betti' num2str(betti-1) ' Band Position=' bandname2{bp}]); clim([0 1]);
ax = gca;
ax.XTick = 1:size(ttestmats, 2);
ax.YTick = 1:size(ttestmats,2);
xticklabels(verbessStates);
yticklabels(verbessStates);
end
saveas(gcf, sprintf('%s9ratTTest_bp=%.2f_betti=%.2f.png', figurepath, bp, betti))
end
end
