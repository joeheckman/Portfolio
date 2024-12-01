---CODE SAMPLE. SCRIPT USED TO PROCESS RAW EEG FOR TOPOLOGICAL DATA ANALYSIS. 'WPLI (WEIGHTED PHASE LAG INDEX)' WAS THE RELEVANT METRIC EXTRACTED FROM THE RAW SIGNAL THAT THE TOPOLOGICAL DATA ANALYSIS WAS PERFORMED ON.---







%%THIS IS THE SCRIPT THAT DID THE WPLI FOR ALL OF THE TDA RESULTS AS OF
%%6.21.24



datapath = 'D:\Joe\TDA_for_rat_data_summer_2024\ns2Files_10_rats\';
savepath = 'D:\Joe\TDA_for_rat_data_summer_2024\EEG_data\';

files = ls(fullfile(datapath, '*.ns2'));

for file = 1:size(files,1)
filename = strtrim(files(file,:));
openNSx('report','read',[datapath files(file,:)],'uV', 'p:double');
EEG = NS2.Data;
save([savepath 'signal' strtrim(filename(1:end-4)) '.mat'], 'EEG', '-v7.3');
end

clc; clear;

%%calculate connectivity, probably wpli

%%%%%% load files %%%%%%%%%%%%%%%%%%%%%%
EEG_path = 'D:\Joe\TDA_for_rat_data_summer_2024\EEG_data';
savepath = ('D:\Joe\TDA_for_rat_data_summer_2024\PyCliqueTop_2023-main\wpli_untrimmed\');
addpath(EEG_path);
addpath(genpath(savepath));

rawEEGfiles = dir(EEG_path);
rawEEGfiles = rawEEGfiles(3:end);


bandname={'d', 't', 'a', 'b', 'g1', 'g2', 'g3', 'b'};
bandname2={'DELTA', 'THETA', 'ALPHA', 'BETA', 'GAMMA1', 'GAMMA2','GAMMA3', 'BROAD'};

band_freq = [1,4;4,10;10,15;15,25;25,55;85,125;125,155;1,155]; % [0.1 2:2:26];

sf_new = 500;  % desired sampling frequency (Hz) post downsampling
old_sf = 1000;



for n = 8:size(rawEEGfiles,1)
load(rawEEGfiles(n).name)
disp(['working with' rawEEGfiles(n).name])

%%this is a cut from some previous analysis that sets up EEGlab for raw
%%time series analysis, nice
ch = size(EEG,1)-2; %number of channels

EEG = double(EEG); 
EEG = resample(EEG',sf_new,old_sf)'; % down-sampling
EEG_cut = EEG(1:30,:); % EEG segments we're using, removing reference and ground
EEG_ref = EEG(31,:); % EEG reference channel
EEG_ref_sig = (EEG_cut-EEG_ref) - sum(EEG_cut-EEG_ref)/ch; % make average reference signals


%%%%%% set variables to use eeglab functions %%%%%%%%%%%
EEG2 = struct();
EEG2.data = EEG_ref_sig;
EEG2.srate = sf_new;
EEG2.chanlocs = []; %%%hmm
EEG2.nbchan = size(EEG_cut,1);
EEG2.trials = 1;
EEG2.event = [];
EEG2.pnts = size(EEG_cut,2);
EEG2.etc = [];

winsize = (20)*sf_new; %this is the size of the data you are analyzing
winmove = (10)*sf_new; %this is the amount the moving window moves for each analysis, NOT the amount of overlap
wmax = floor((size(EEG2.data',1)-winsize+winmove)/(winmove));



scalingFactor = 3600*sf_new/(winmove) ; %(samplingrate*seconds)/(samplingrate*seconds) => nondeminsionalized scaling factor
analysis_windows = {.4,.5;1.4,1.5; 2.35,2.45;3.1,3.2};%%in hours
analysis_states = {'N_W', 'S_1', 'S_2', 'RoC'};
analysis_frames = cellfun(@(x) x*scalingFactor,  analysis_windows);



WPLI = cell(size(band_freq,1),1);


for bp = 1:size(band_freq,1)
filtEEG = pop_eegfiltnew(EEG2, band_freq(bp,1),band_freq(bp,2), 3300); % filtering
h_data = filtEEG.data';
htransformdata = hilbert(h_data);

for row = 1:size(analysis_frames,1)
analyse_these = analysis_frames(row,1):1:analysis_frames(row,2);
analyse_these = round(analyse_these);
disp(['generating wpli for state ' analysis_states{row}]);

for w = analyse_these(1):analyse_these(end)
disp(['n=' num2str(n) ',bp=' num2str(bp) ',w=' num2str(w)]);
a_sig = htransformdata(winmove*(w-1)+1:winmove*(w-1)+winsize,:); % window-segmented data to be analyzed. 
Wpli = w_PhaseLagIndex2(a_sig); % calculate wpli
wpli(:,:,mod(w,analyse_these(1))+1) = Wpli; %%%CAREFUL WITH THISLINE

end
save(sprintf(['%swpli_bp=%.2f_' analysis_states{row} rawEEGfiles(n).name ],savepath, bp), 'wpli' )

end
end
end