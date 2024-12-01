%%%CODE SAMPLE. THIS SCRIPT WAS USED TO GENERATE FIGURES OF VARIANCE IN ORDER PARAMETER CURVES FOR KURAMOTO MODEL SIMULATIONS AS IT SYNCHRONIZED AND DESYNCHRONIZED. 
%%%FORWARD AND BACKWARD CAN BE COMPARED WITH THIS, LAST TEST DID NOT INCLUDE BACKWARD


%%last use: FORWARDS TEST ON HIGHER RES THROUGH SYNCH. 

clc; clear;

%%set datapath to time series
Datapath= 'C:\Kura_Or_Orstd_pcf_curve_shape_test_higher_res_take_1\';
addpath(Datapath)
filedir = Datapath;

maxindexsavepath = 'D:\Joe\Zauberbaum\DataGeneration\BetterKuraResults\Kura_maxindex_higher_res_take_1\';
% figure_savepath = 'D:\Joe\Zauberbaum\DataGeneration\BetterKuraResults\reduced_net_figures_1.26\';

%%number of iterations, coupling strength vector, target average frequency gap
%%calculated by 2013 paper we are following
steps = 2000;
cores = 16; sets = 5; K_vec = [linspace(0,.25,steps) linspace(.25,0,steps)]; freq_gap = .35:.01:.44;
Filelist = dir(filedir);
Filelist = Filelist(3:numel(Filelist));

files = ls(fullfile(filedir, '*.mat'));
fileCellArray = cellstr(files);


for FG= 1:length(freq_gap)
    searchText = sprintf('freq_gap=%.2f', freq_gap(FG));

    %%find all matching files for forward iteration order param curve
    OR = [];
    ORSTD = [];
    for a = 1:steps
        disp(a)
        tic
        parfor b = 1:(cores * sets)

            searchText2 = sprintf('state=%.2f', a);
            searchText3 = sprintf('iteration=%.2f',b);
            matchingFile = fileCellArray(contains(fileCellArray, searchText) ...
                & contains(fileCellArray, searchText2) ...
                & contains(fileCellArray, searchText3));
%             disp('Matching Files:');
            x = matchingFile{1};
            y = load(x);
            OR = [OR;y.Or_vec];
            ORSTD = [ORSTD;y.Or_std_vec];
            % Display the matching files
%             disp(matchingFile);
        end
        toc
    end
    Or_data_reshaped = reshape(OR,cores * sets,steps); %%reshapes to matrix col = iterations, row= index of K_vec
    Or_std_data_reshaped = reshape(ORSTD,cores * sets,steps);%%same with orstd
    save(sprintf('highres_kura_test_Fg=%.3f.mat', freq_gap(FG)), 'Or_std_data_reshaped');
end

 


%%same thing backwards
    ORbackward = [];
    ORSTDbackward =[];
    for a = steps+1:2*steps
        parfor b = 1:(cores * sets)
            searchText2 = sprintf('state=%.2f', a);
            searchText3 = sprintf('iter_no=%.2f',b);
            matchingFile = fileCellArray(contains(fileCellArray, searchText) ...
                & contains(fileCellArray, searchText2)...
                & contains(fileCellArray, searchText3));
            disp('Matching Files:');
            x = matchingFile{1};
            y = load(x);
            ORbackward = [ORbackward;y.Or_vec];
            ORSTDbackward = [ORSTDbackward;y.Or_std_vec];

            % Display the matching files
            disp(matchingFile);
        end
    end
    Or_reverse_data_reshaped = reshape(ORbackward,cores * sets,steps); %%reshapes to matrix col = iterations, row= index of K_vec
    Or_reverse_std_data_reshaped = reshape(ORSTDbackward,cores * sets,steps);%%same with orstd
  


    %%find the maximum index of OR_STD as critical point of syncrhonization
    for m = 1:size(Or_std_data_reshaped,1)
        [maxValue(m), Index(m)] =max(Or_std_data_reshaped(m,:));
%         [reversemaxValue(m), reverseIndex(m)] = max(Or_reverse_std_data_reshaped(m,:));
        disp(['maximum orstd value at ' num2str(Index(m)) ' value is ' ...
            num2str(max(Or_std_data_reshaped(m,:)))]);

%         disp(['reverse, maximum orstd value at ' num2str(reverseIndex(m)) ' value is ' ...
%             num2str(max(Or_reverse_std_data_reshaped(m,:)))]);
    end
    % save(sprintf('%sMaxPCFKvecindex_freqgap=%.3f.mat', maxindexsavepath, freq_gap(FG)), 'Index');
    save(sprintf('%sMaxPCFKvecindex_freqgap=%.3f.mat', maxindexsavepath, freq_gap(FG)), 'Index', 'reverseIndex'); %%%saves index in order of iter_no


    %%align curves for figure
    shift = Index-max(Index);
%     reverseShift = reverseIndex-max(reverseIndex);
    Extended_mat = [Or_std_data_reshaped, zeros(size(Or_std_data_reshaped,1), abs(min(shift)))];
%     reverseExtended_mat = [Or_reverse_std_data_reshaped, zeros(size(Or_reverse_std_data_reshaped,1), abs(min(reverseShift)))];
    for o = 1:size(Or_std_data_reshaped,1)
        Extended_mat(o,:) = circshift(Extended_mat(o,:),abs(shift(o)));
%         reverseExtended_mat(o,:) = circshift(reverseExtended_mat(o,:),abs(reverseShift(o)));
        [maxValueShift(o), Index2(o)] = max(Extended_mat(o,:));
%         [reversemaxValueShift(o), reverseIndex2(o)] = max(reverseExtended_mat(o,:));

        disp(['maximum aligned orstd value at ' num2str(Index2(o))])
%         disp(['maximum aligned reverse orstd value at ' num2str(reverseIndex2(o))])

    end

    %%aligned curves for forward
    figure();
    for figure_num = 1:cores

        ylabel('Variance in Order Parameter'); ylim([0 0.4]);
        xlim([abs(min(shift)) steps+abs(min(shift))]);
        ylim([0 .275])
        set(gca,'Xticklabel', [])
        set(gca,'XTick', [])
        title(['Variance in order parameter aligned to max values. Frequency gap is ' num2str(freq_gap(FG))])
        set(gca,'fontsize',8, 'FontName', "Times New Roman")
        plot(Extended_mat(figure_num,:),'marker', 'none', 'linestyle', '-', 'linewidth',.6, 'color', [0 0 0])
        hold on
    end


    %%aligned curves for backward
    figure();
    for figure_num = 1:iter_no
        plot((reverseExtended_mat(figure_num,end:-1:1)),'marker', 'none', 'linestyle', '-', 'linewidth',.6, 'color', [0 0 0])
        ylabel('Variance in Order Parameter'); ylim([0 0.4]);
        xlim([abs(min(reverseShift)) steps+abs(min(reverseShift))]);
        ylim([0 .275])
        set(gca,'Xticklabel', [])
        set(gca,'XTick', [])
        title(['Backwards, variance in order parameter aligned to max values. Frequency gap is ' num2str(freq_gap(FG))])
        set(gca,'fontsize',8, 'FontName', "Times New Roman")
        hold on

        newFigureSize = [0, 0, 1000, 1300];
        set(gcf, 'Position', newFigureSize);

    end


end
iter_no = 16;

tiledlayout(4,4);
tiledlayoutPosition = get(gcf, 'Position');
tiledlayoutPosition(3) = tiledlayoutPosition(3) * 1.5; % Increase width
tiledlayoutPosition(4) = tiledlayoutPosition(4) * 1.5; % Increase height
set(gcf, 'Position', tiledlayoutPosition);
for or_fig2 = 1:iter_no
    nexttile
    hold on

    plot((Or_data_reshaped(or_fig2,:)),'marker', 'none', 'linestyle', '-', 'linewidth',.6, 'color', [0 0 0])
%     plot((Or_reverse_data_reshaped(or_fig2,end:-1:1)),'marker', 'none', 'linestyle', '-', 'linewidth',.6, 'color', [1 0 1])
    ylabel('OP'); ylim([0 1])
    set(gca,'fontsize',6, 'FontName', "Times New Roman")
    set(gca,'fontsize',6, 'FontName', "Times New Roman")
    xlim([1 steps]);
    set(gca,'Xticklabel', [])
    set(gca,'XTick', [])
    title(['Order Param curves. Frequency gap is ' num2str(freq_gap(FG))])
    set(gca,'fontsize',6, 'FontName', "Times New Roman")


end
saveas(gcf,sprintf('%sOr_orstd_curves_with_max_index_freqGap=%.3f.png',figure_savepath, freq_gap(FG)))
end

close all;




