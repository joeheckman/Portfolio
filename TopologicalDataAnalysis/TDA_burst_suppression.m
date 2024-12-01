%%%CODE SAMPLE. ANOTHER SCRIPT USED TO CREATE VIDEO SHOWING TOPOLOGICAL INFORMATION OF RAT EEG WHILE UNDER DEEP ANESTHESIA




%%%see what's happening with burst suppresion. This script and one other
%%%are used for getting Wpli 
load('R081021A_TTX_sevo_210907.mat')
 figure(); clf
  spectrogram(EEG(1,:),fs*.8,[],0:.05:35,fs,  "yaxis"); 
 numTicks = 200; 
  xLimits = xlim;
  xTicks =linspace(xLimits(1), xLimits(2), numTicks); 
  xticks(xTicks);
  xticklabels(xticks)
    colormap jet; 
    clim([0 35]); 
  


%%processing. Run eeglab
load('R081021A_TTX_sevo_210907.mat'); fs = 1000; %%careful with fs. It was given in data file
addpath('C:\Program Files\MATLAB\packages\eeglab2024.0\')
savepath = ('D:\Joe\TDA_for_rat_data_summer_2024\PyCliqueTop_2023-main\Data_Network\');
bandname={'d', 't', 'a', 'b', 'g1', 'g2', 'g3', 'b'};
bandname2={'DELTA', 'THETA', 'ALPHA', 'BETA', 'GAMMA1', 'GAMMA2','GAMMA3', 'BROAD'};

band_freq = [1,4;4,10;10,15;15,25;25,55;85,125;125,155;1,155]; % [0.1 2:2:26];


%%this is a cut from some previous analysis that sets up EEGlab for raw
%%time series analysis, nice
ch = size(EEG,1)-2; %number of channels

EEG = double(EEG); 
% EEG = resample(EEG',sf_new,old_sf)'; maybe I won't downsample this time
EEG_cut = EEG(1:30,:); % EEG segments we're using, removing reference and ground
EEG_ref = EEG(31,:); % EEG reference channel
EEG_ref_sig = (EEG_cut-EEG_ref) - sum(EEG_cut-EEG_ref)/ch; % make average reference signals


%%%%%% set variables to use eeglab functions %%%%%%%%%%%
EEG2 = struct();
EEG2.data = EEG_ref_sig;
EEG2.srate = fs;
EEG2.chanlocs = []; %%%hmm
EEG2.nbchan = size(EEG_cut,1);
EEG2.trials = 1;
EEG2.event = [];
EEG2.pnts = size(EEG_cut,2);
EEG2.etc = [];




win = 3
winsize = (win)*fs; %this is the size of the data you are analyzing
winmove = (1/200)*fs; %this is the amount the moving window moves for each analysis, NOT the amount of overlap


analysis_windows = {.48,52};%%in hours. Visual inspection
analysis_frames = cellfun(@(x) x*3600*fs,  analysis_windows);

%%the generated result here goes to Data_network. Run Spyder script, 
% It ends up in TDA_results.
bp = 2 %%theta 
filtEEG = pop_eegfiltnew(EEG2, band_freq(bp,1),band_freq(bp,2), 3300); % filtering
h_data = filtEEG.data';
htransformdata = hilbert(h_data);

wpli = [];
for w = analysis_frames(1):winmove:analysis_frames(2)
disp(['time=' num2str(w/(3600*fs))])
a_sig = htransformdata(w:w+winsize,:); % window-segmented data to be analyzed. 
Wpli = w_PhaseLagIndex2(a_sig); % calculate wpli
wpli = cat(3,wpli,Wpli);
end
save(sprintf(['%swpli_bp=%.2f_' 'winsize=%.2f_overlap=%.4f_highres_BurstSuppression.mat' ],savepath, bp, win, winmove/fs), 'wpli' )
% clear wpli



%%visualize
win = 3; fs = 1000;
winsize = (1/200)*fs; %this is the size of the data you are analyzing
winmove = (20)*fs; %t
datapath = 'D:\Joe\TDA_for_rat_data_summer_2024\PyCliqueTop_2023-main\TDA_results';
addpath(datapath);
savepath2 = 'D:\Joe\TDA_for_rat_data_summer_2024\highTimeRes_TDA\videotest\';
files = ls(fullfile(datapath, '*.mat'));
files = cellstr(files);

ylims={50,50,50,50,50,30};

    winsize_searchtext = sprintf('winsize=%.2f', win);
    overlap_searchtext = sprintf('overlap=%.2f', winmove/fs);
    file = files(contains(files, winsize_searchtext) & contains(files, overlap_searchtext));
load(file{1,1}) 



% v = VideoWriter([savepath2 'throughaLOCtest'], 'MPEG-4');
v.FrameRate = 5;
open(v);
new_position = [0,0,1200,1200];
        set(figure, 'Position', new_position); 
       for t = 1:size(BettiCurves,1)
       clf
    
            for i = 1:(size(BettiCurves,3))
                subplot(size(BettiCurves,3),1,i); ylabel('# Cycles'); xlabel(['BN' num2str(i-1)])
        % if i == 1
        %     % title(['time elapsed=' num2str((t-1)/v.FrameRate) ' seconds']);
        % end
        %         ylim([0 ylims{i}]);
                    hold on
                    plot(edgeDensity, BettiCurves(t,:,i), 'Color', [1 0 1])
                    set(gca, 'XTickLabel', get(gca, 'XTickLabel'), 'FontSize', 11)
           
            end
              figHandle = gcf;
              set(findall(figHandle, '-property', 'FontName'), 'FontName', 'Times New Roman');  
            drawnow
            hold on
 
frame = getframe(gcf);
writeVideo(v,frame);
       end
