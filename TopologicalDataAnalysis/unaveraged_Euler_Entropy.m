---CODE SAMPLE. THIS SCRIPT GENERATES 'EULER ENTROPIES' FOR TOPOLOGICAL SIGNATURE IN RAT EEG. EULER ENTROPY IS USED FOR IDENTIFYING PHASE TRANSITIONS IN NET WORK TOPOLOGY. THIS SCRIPT IS ADJACENT TO FIGURES INCLUDED IN THE PORTFOLIO BUT NONETHELESS MEANINGFUL.---



bettispath = 'D:\Joe\TDA_for_rat_data_summer_2024\PyCliqueTop_2023-main\TDA_results_9rat\';
figurepath = 'D:\Joe\TDA_for_rat_data_summer_2024\PyCliqueTop_2023-main\EE_result_9rat\';
addpath(bettispath);
files = ls(fullfile(bettispath, '*.mat'));
files = cellstr(files);
bandname2={'DELTA', 'THETA', 'ALPHA', 'BETA', 'GAMMA1', 'GAMMA2','GAMMA3', 'BROAD'};
band_freq = {'1-4hz','4-10hz','10-15hz','15-25hz','25-55hz','85-125hz','125-155hz','1-155hz'};
states = {'N_W','RoC', 'S_1', 'S_2'};
verbessStates = {'Normal Wake', 'Recovery', 'Sevo Epoch 1', 'Sevo Epoch 2'};
state_colors = {'Red', 'magenta', 'Blue', 'black'};
ylims={30,25,20,18,12,4, 1234523613613};


for bp = 2
    bp_searchtext = sprintf('%.2f', bp);
    betFiles = files(contains(files,bp_searchtext));
    betFiles = reshape(betFiles, numel(betFiles)/numel(states), []);
    for p = 1:size(betFiles,2)
        for m = 1:size(betFiles,1)
    figure();
            load(betFiles{m,p})
            disp(betFiles{m,p})
            avgBettiCurves = mean(BettiCurves,1);
            sqBettiCurves = squeeze(avgBettiCurves);
            EE = Betti2EulerEntropy(sqBettiCurves);
    plot(edgeDensity,EE, 'color', [1 0 1]);title(['bp=' num2str(bp_searchtext) ' ' verbessStates{p}]);
    ylabel('Euler Entropy'); xlabel('Edge Density'); ylim([-5 5])
    set(gca,'FontName', 'Times New Roman')
    saveas(gcf, [sprintf('%s9ratUnaveragedEE_bp=%.2f_', figurepath, bp) 'rat#=' num2str(m) states{p} '.png'])

        end
    end
end