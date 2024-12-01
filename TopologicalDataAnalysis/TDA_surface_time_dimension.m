---CODE SAMPLE. THIS SCRIPT VISUALIZES TOPOLOGICAL PATTERNS IN RAT EEG THROUGH TIME WHERE TIME ELAPSED IS A DIMENSION OF THE SURFACE (SEE RESULT IN PORTFOLIO). HIGHLY EFFECTIVE.---





%%%make a surface of High Time Res TDA with time axis as a dimension, end
%%%goal being for each betti number for each state. 
datapath='D:\Joe\TDA_for_rat_data_summer_2024\PyCliqueTop_2023-main\TDA_results_highTimeRes_4epochs_9rat\';
addpath(datapath)
%IMPORTANT: WINSIZE FOR ABOVE DATAPATH: 4. OVERLAP: 1/50. WASN'T DOCUMENTED
%IN FILENAMES. 


% load('TDA_bp=2.00_winsize=3.00_overlap=0.0050_highres_BurstSuppression.mat')
% (*for burst suppression tests*)


result_path='D:\Joe\TDA_for_rat_data_summer_2024\highTimeRes_TDA\9rat_TDA_surfaces\';
analysis_states = {'NormalWake', 'Recovery', 'Sevo1', 'Sevo2'};
files = cellstr(ls(fullfile(datapath, '*mat')));
files = reshape(files,[],size(analysis_states,2));
winsize = 4; 
overlap = 1/50; 
zlims={100,100,100,100,100,100};

 
%%X = edgeDensity;
%%Y = time dimension
for i = 1:size(files,2)
    for p = 4:size(files,1) %%just look at 3 rats
       
    load(files{p,i})
    disp(files{p,i})
Y = 1:size(BettiCurves,1);

        for m = 1:size(BettiCurves,3)
            new_position = [0.2,386,2044,475];
            set(figure, 'Position', new_position);
surface = surf(edgeDensity,Y,BettiCurves(Y(1):Y(end),1:end,m));
surface.EdgeColor = 'none'; surface.FaceColor = 'texturemap'; colormap jet
view(-53,51); zlabel('# cycles'); 


zlim([0 100]); 
zticks(0:10:200);

if i == 2
    zlim([0 200]);
    zticks(0:20:200);
end

xlabel('Edge Density'); 
ylabel('Time Elapsed');yticks(1:300:max(Y)); ylim([1 max(Y)]); yticklabels((yticks-Y(1))*overlap);
title(['Rat#=' num2str(p) ' ' analysis_states{i} ', BettiNumber=' num2str(m-1)])
set(gca,'FontName', 'Times New Roman')
saveas(gcf, [sprintf('%sTDAsurface_Rat%.0f_Betti=%.0f_state=', result_path, p, m) analysis_states{i} '.fig'])
saveas(gcf, [sprintf('%sTDAsurface_Rat%.0f_Betti=%.0f_state=', result_path, p, m) analysis_states{i} '.png'])
        end
        close all
    end
end