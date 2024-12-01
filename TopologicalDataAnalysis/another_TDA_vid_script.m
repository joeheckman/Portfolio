---CODE SAMPLE. THIS IS A SCRIPT USED TO CREATE VIDEOS SHOWING THE EVOLUTION OF TOPOLOGICAL INFORMATION IN RAT EEG AS ITS CONSCIOUS STATE IS BEING PERTURBED.---







%%use recovery of consciuosness period to see topological changes. 
%%spectrogram
eeg_path = 'D:\Joe\TDA_for_rat_data_summer_2024\highTimeRes_TDA\eegsig';
addpath(eeg_path);
fs = 1000;
sigdata = dir(eeg_path);
sigdata=sigdata(3:end);
for m = 6
load(sigdata(m).name)
disp(sigdata(m).name)
figure();
  spectrogram(EEG(1,:),fs*.8,[],0:.05:35, fs,  "yaxis"); 
  clim([0 35])
  colormap jet
  title(sigdata(m).name)
  numTicks = 500; 
  xLimits = xlim;
  xTicks =linspace(xLimits(1), xLimits(2), numTicks); 
  xticks(xTicks);
  % saveas(gcf, sprintf(['%s' sigdata(m).name(1:end-5) 'spectro' '.png'], specpath))
end


%%processing. Run eeglab
eegpath = 'D:\Joe\TDA_for_rat_data_summer_2024\EEG_data\eegsig\';
addpath(eegpath)
addpath('C:\Program Files\MATLAB\packages\eeglab2024.0\')
savepath = ('D:\Joe\TDA_for_rat_data_summer_2024\PyCliqueTop_2023-main\Data_Network\');
bandname={'d', 't', 'a', 'b', 'g1', 'g2', 'g3', 'b'};
bandname2={'DELTA', 'THETA', 'ALPHA', 'BETA', 'GAMMA1', 'GAMMA2','GAMMA3', 'BROAD'};

band_freq = [1,4;4,10;10,15;15,25;25,55;85,125;125,155;1,155]; % [0.1 2:2:26];

rawEEGfiles = ls(fullfile(eegpath, '*.mat'));
rawEEGfiles = cellstr(rawEEGfiles);


for n = 1:size(rawEEGfiles,1)
load(rawEEGfiles{n,1})
disp(['working with ' rawEEGfiles{n,1}])

%%this is a cut from some previous analysis that sets up EEGlab for raw
%%time series analysis, nice

fs = 1000;
ch = size(EEG,1)-2; %number of channels
EEG = double(EEG); 
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


win = 4
winsize = (win)*fs; %this is the size of the data you are analyzing
winmove = (1/50)*fs; %this is the amount the moving window moves for each analysis, NOT the amount of overlap


analysis_states = {'Normal', 'Sevo1', 'Sevo2', 'RoC'};
analysis_windows = {.4,.425;1.4,1.425; 2.35,2.375;3.1,3.125};%%in hours
analysis_frames = cellfun(@(x) x*3600*fs,  analysis_windows);

bp = 2
filtEEG = pop_eegfiltnew(EEG2, band_freq(bp,1),band_freq(bp,2), 3300); % filtering
h_data = filtEEG.data';
htransformdata = hilbert(h_data);

for row = 1:size(analysis_frames,1)
wpli = [];
for w = analysis_frames(row,1):winmove:analysis_frames(row,2)
disp(['time=' num2str(w/(3600*fs))])
a_sig = htransformdata(w:w+winsize,:); % window-segmented data to be analyzed. 
Wpli = w_PhaseLagIndex2(a_sig); % calculate wpli
wpli = cat(3,wpli,Wpli);
end
save(sprintf(['%sHighTimeResWpli_bp=%.2f_' analysis_states{row} rawEEGfiles{n,1} ],savepath, bp), 'wpli' )
end
end



%%visualize
win = 4; fs = 1000;
winsize = (win)*fs; %this is the size of the data you are analyzing
winmove = (1/50)*fs; %t
datapath = 'D:\Joe\TDA_for_rat_data_summer_2024\PyCliqueTop_2023-main\TDA_results';
addpath(datapath);
savepath2 = 'D:\Joe\TDA_for_rat_data_summer_2024\highTimeRes_TDA\videotest\';
files = ls(fullfile(datapath, '*.mat'));
files = cellstr(files);

ylims={50,50,50,50,50,30};

   searchtext=('TimeRes');
    files = files(contains(files, searchtext));

    for m = 1:size(files,1)
        load(files{m,1})

        v = VideoWriter([savepath2 files{m,1}(24:end) 'asdfasdf'], 'MPEG-4');
        v.FrameRate = 50;
        open(v);
        new_position = [0,0,1200,1200];
        set(figure, 'Position', new_position);
        for t = 1:size(BettiCurves,1)
            clf

            for i = 1:size(BettiCurves,3)
                subplot(size(BettiCurves,3),1,i); ylabel('# Cycles'); xlabel(['BN' num2str(i-1)])
                if i == 1
                    title(['time elapsed=' num2str((t-1)/v.FrameRate) ' seconds']);
                end
                ylim([0 ylims{i}]);
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

    end