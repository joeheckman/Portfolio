---CODE SAMPLE. SCRIPT COMPARES A SET OF TOPOLOGICAL METRICS (BETTI NUMBERS 1-5 AND EULER ENTROPY) BETWEEN ALTERED AND UNALTERED STATES OF CONSCIOUSNESS IN RAT EEG.---


%%figure for different conscious states overlayed
bettispath = 'D:\Joe\TDA_for_rat_data_summer_2024\PyCliqueTop_2023-main\TDA_results_9rat\';
figurepath = 'D:\Joe\TDA_for_rat_data_summer_2024\PyCliqueTop_2023-main\TDA_visual_result\';
figurepath2 = 'D:\Joe\TDA_for_rat_data_summer_2024\PyCliqueTop_2023-main\test_result\';
addpath(bettispath);


files = ls(fullfile(bettispath, '*.mat'));
files = cellstr(files);
bandname2={'DELTA', 'THETA', 'ALPHA', 'BETA', 'GAMMA1', 'GAMMA2','GAMMA3', 'BROAD'};
band_freq = {'1-4hz','4-10hz','10-15hz','15-25hz','25-55hz','85-125hz','125-155hz','1-155hz'};
states = {'N_W', 'S_1', 'S_2', 'RoC'};
verbessStates = {'Normal Wake', 'Sevo Epoch 1', 'Sevo Epoch 2', 'Recovery'};
state_colors = {'Red', 'magenta', 'Blue', 'black'};
ylims={30,25,20,18,12,4, 1234523613613};


for bp = 1:8
    bp_searchtext = sprintf('%.2f', bp);
    betFiles = files(contains(files,bp_searchtext));
    betFiles = reshape(betFiles, numel(betFiles)/numel(states), []);

    if bp ==1 %%find min edge density for interp1

        findMaxED = [];

        for m = 1:size(betFiles,1)


            load(betFiles{m,1})
            disp(betFiles{m,1})
            avgBettiCurves = mean(BettiCurves,1);
            sqBettiCurves = squeeze(avgBettiCurves);
            EE = Betti2EulerEntropy(sqBettiCurves);
            findMaxED = [findMaxED; numel(EE)];

        end
        maxEdgeDensity = linspace(0,1,max(findMaxED));
    end



    %%with max edge density, create common edge density values (because of
    %%different number of missing channels and network size determines edge
    %%density) and make first figure.

    for m = 1:size(betFiles,1)
        new_position = [0,0,1200,1200];
        set(figure, 'Position', new_position);
        for p = 1:size(betFiles,2)



            load(betFiles{m,p})
            disp(betFiles{m,p})
            avgBettiCurves = mean(BettiCurves,1);
            sqBettiCurves = squeeze(avgBettiCurves);
            EE = Betti2EulerEntropy(sqBettiCurves);
            interpBettiCurves = interp1(edgeDensity, sqBettiCurves, maxEdgeDensity);


            for i = 1:(size(sqBettiCurves,2)+1)
                subplot(size(BettiCurves,3)+1,1,i); ylabel('# Cycles');
                ylim([0 ylims{i}]);


                if i == 1
                    title([bandname2{bp} ' ' band_freq{bp} ' rat#=' num2str(m)])
                end
                if i < (size(sqBettiCurves,2)+1)

                    hold on
                    % plot(edgeDensity, sqBettiCurves(:,i), 'Color', state_colors{p})
                    plot(maxEdgeDensity, interpBettiCurves(:,i), 'Color', state_colors{p})
                    set(gca, 'XTickLabel', get(gca, 'XTickLabel'), 'FontSize', 11)
                else
                    plot(edgeDensity,EE, 'Color', state_colors{p})
                    title('EE vs. Edge Density'); xlabel('Edge Density', FontSize=15); ylabel('EE');
                    set(gca, 'XTickLabel', get(gca, 'XTicklabel'), 'FontSize', 11)
                    ylim([-5 5]);
                end
            end
            hold on


        end
        figHandle = gcf;
        set(findall(figHandle, '-property', 'FontName'), 'FontName', 'Times New Roman');
        legend(verbessStates)
        saveas(gcf, sprintf(['%s' betFiles{m,1}(1:11) betFiles{m,1}(21:end-12) '.png'], test_figurepath))
    end
end






%%a different figure that gives averaged results across datasets (I don't
%%trust this result at all lol)
bettispath = 'D:\Joe\TDA_for_rat_data_summer_2024\PyCliqueTop_2023-main\TDA_results\';
figurepath = 'D:\Joe\TDA_for_rat_data_summer_2024\PyCliqueTop_2023-main\TDA_visual_result\';
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
    end %found max edge density resolution that everything needs to be interpolated to

    bettisAcrossSets= zeros(numel(maxEdgeDensity), size(sqBettiCurves,2),size(betFiles,1)); %%dataset x edgedensity x betti number
    for m = 1:size(states,2) %%for each state
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
        avgEEAcrossSets(:,m) = mean(interpEE,2);
        avgBettisAcrossSets(:,:,m) = mean(bettisAcrossSets,3); %%avg figure for bettis 0-5 of all sets, third dimension is across states for the figure
    end


 %%%%figure for each state
        new_position = [0,0,1200,1200];   %make a figure for each state
        set(figure, 'Position', new_position);
        for r = 1:size(avgBettisAcrossSets,2)+1 %%entries at each Betti

            subplot(size(avgBettisAcrossSets,2)+1,1,r); ylabel('# Cycles');


            for fig2 = 1:size(avgBettisAcrossSets,3) %%who cares
                if r < size(avgBettisAcrossSets,2)+1
                plot(maxEdgeDensity, avgBettisAcrossSets(:,r,fig2), 'Color', state_colors{fig2})
                 ylabel('# Cycles')
                xlabel(sprintf('Edge Density, BN=%.0f', r-1))
                set(gca, 'XTickLabel', get(gca, 'XTickLabel'), 'FontSize', 11)
                else
                    plot(maxEdgeDensity,avgEEAcrossSets(:,fig2), 'Color', state_colors{fig2});
                    xlabel('Edge Density')
                    ylabel('Euler Entropy')
                end
               
                % ylim([0 ylims{r}]);
                hold on
            end

            if r < size(avgBettisAcrossSets,2)+1
            lim = max(avgBettisAcrossSets(:,r,:));
            ylim([0 max(lim)+1]);
            if r == 1
                title([bandname2{bp} ' ' band_freq{bp} 'All States'])
                legend(verbessStates)
            end
            end
            figHandle = gcf;

            set(findall(figHandle, '-property', 'FontName'), 'FontName', 'Times New Roman');
        end
        saveas(gcf, sprintf(['%s' betFiles{m,1}(1:11) 'avgBettiCurves_allstates' '.png'], figurepath2))
end