---CODE SAMPLE. A FAILED TEST FOR STARTING KURAMOTO MODEL SIMULATIONS AT THE MOMENT BEFORE SYNCHRONIZATION TO LIMIT COMPUTATIONAL TIME AND STOARGE ---

datapath = 'C:\Kura_Orstd_truncated_test\';
sig_datapath = 'D:\Joe\Zauberbaum\DataGeneration\BetterKuraResults\Kura_sig_truncated_test\';
addpath(datapath)
addpath(sig_datapath)

addpath('D:\Joe\Zauberbaum\AdjustFrequencyGapNetworks\pcf_curve_shape_test')
Kura_truncated_sig_savepath = 'D:\Joe\Zauberbaum\DataGeneration\BetterKuraResults\Kura_sig_truncated_result_test\';
Kura_truncated_or_orstd_savepath = 'C:\Kura_Orstd_truncated_result_test\';
filelist = dir(datapath);
filelist = filelist(3:end);
data = ls(fullfile(datapath, '*.mat'));
sig_data = ls(fullfile(sig_datapath, '*.mat'));
Orstddata = cellstr(data);
sigdata = cellstr(sig_data);
target_FG = .35;

steps = 100;
cores = 16; sets = 1; K_vec = linspace(0,.2,steps); freq_gap = .35; 


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
            matchingFile = Orstddata(contains(Orstddata, searchText) ...
                & contains(Orstddata, searchText2) ...
                & contains(Orstddata, searchText3));

            x = matchingFile{1};
            y = load(x);
            OR = [OR;y.Or_vec];
            ORSTD = [ORSTD;y.Or_std_vec];
      
        end
        toc
    end
    Or_data_reshaped = reshape(OR,cores * sets,steps); %%reshapes to matrix col = iterations, row= index of K_vec
    Or_std_data_reshaped = reshape(ORSTD,cores * sets,steps);%%same with orstd
%     save(sprintf('highres_kura_test_Fg=%.3f.mat', freq_gap(FG)), 'Or_std_data_reshaped');
end

 

figure();
for fig = 1:size(Or_data_reshaped,1)
    plot(K_vec, Or_data_reshaped(fig,:))
    hold on
end

%%say Ki =30 is where we want to start 



for i = 1:length(target_FG)
load(sprintf('NaturalFrequencies_and_MAT_FreqGap=%.3f.mat', target_FG(i)));
disp(['running simulation for target FG ' num2str(target_FG(i))])

sets = 1; Ki_start = 30; %% from figures above, start of griffiths phase

T=0:0.01:80;
dt=T(2)-T(1);
transient_time=40; cores = 16; %% 16 cores
iterations = sets * cores;
tic
noise_std=0; %%%no noise in 2013 paper we are following 10.12

or =zeros(length(K_vec),iterations);
or_std = zeros(length(K_vec),iterations);
K_vec_length = length(K_vec);


searchText_Ki_start = sprintf('state=%.2f', Ki_start-1);


for set = 1:sets
parfor iter = 1 + ((set-1)*cores): set * cores
    searchText_iter = sprintf('iteration=%.2f',iter);
    Ki_start_file = sigdata(contains(sigdata, searchText_Ki_start)...
        & contains(sigdata, searchText_iter));
    x = Ki_start_file{1};
    Ki_start_sig = load(x);
    w = W_save(iter,:);
    MAT0 = MAT(:,:,iter);
    N=size(MAT,1);
    sig=zeros(N,length(T)-transient_time*1/dt);
% for Ki=1:K_vec_length
 for Ki = Ki_start:K_vec_length %%%first ripser test 3.29
    k=K_vec(Ki);
    MAT2=k*MAT0;
    % noise
    noisein=noise_std*randn(N,length(T));
    if Ki ==Ki_start
%    start = 2*pi*rand(N,1)-pi;
     start = Ki_start_sig.sig(:,end);
   
    else 
   start = phi(:,end);
    end

    % simulation
    [phi] = Kuramoto_noisein(MAT2,T,w',noisein,start);
    phi=phi(:,transient_time*1/dt+1:end);
    [or(Ki,iter),op,or_std(Ki,iter),or_t] = OrderParameter2(sin(phi)');
    Or_vec=or(Ki,iter);
    Or_std_vec=or_std(Ki,iter);

    sig(:,:)=phi;
    
    t=T(1:5:end-transient_time*1/dt);
    sig=sin(sig(:,:)); 
    fname=sprintf('%sKura_OR_ORSTD_freq_gap=%.4f_k=%.6f_iteration=%.3fstate=%.2f.mat', ...
        Kura_truncated_or_orstd_savepath, target_FG(i), k, iter, Ki);
     fname2=sprintf('%sKura_Sig_freq_gap=%.4f_k=%.6f_iteration=%.3fstate=%.2f.mat', ...
         Kura_truncated_sig_savepath,target_FG(i),k,iter,Ki)
     parsave_kura_data(fname,t,w,k,MAT2,Or_std_vec,Or_vec)
     parsave_kura_data_signal(fname2,sig)

end
end
toc
end
end


clc;
clear;



datapath = 'C:\Kura_Orstd_truncated_result_test\';
sig_datapath = 'D:\Joe\Zauberbaum\DataGeneration\BetterKuraResults\Kura_sig_truncated_result_test\';
addpath(datapath)
addpath(sig_datapath)

filelist = dir(datapath);
filelist = filelist(3:end);
data = ls(fullfile(datapath, '*.mat'));
sig_data = ls(fullfile(sig_datapath, '*.mat'));
Orstddata = cellstr(data);
sigdata = cellstr(sig_data);
target_FG = .35; Ki_start = 30;

steps = 100;
cores = 16; sets = 1; K_vec = linspace(0,.2,steps); freq_gap = .35; 

for FG= 1:length(freq_gap)
    searchText = sprintf('freq_gap=%.2f', freq_gap(FG));

    %%find all matching files for forward iteration order param curve
    OR = [];
    ORSTD = [];
    for a = Ki_start:steps
        disp(a)
        tic
        parfor b = 1:(cores * sets)

            searchText2 = sprintf('state=%.2f', a);
            searchText3 = sprintf('iteration=%.2f',b);
            matchingFile = Orstddata(contains(Orstddata, searchText) ...
                & contains(Orstddata, searchText2) ...
                & contains(Orstddata, searchText3));

            x = matchingFile{1};
            y = load(x);
            OR = [OR;y.Or_vec];
            ORSTD = [ORSTD;y.Or_std_vec];
      
        end
        toc
    end
    Or_data_reshaped = reshape(OR,cores * sets,steps-Ki_start+1); %%reshapes to matrix col = iterations, row= index of K_vec
    Or_std_data_reshaped = reshape(ORSTD,cores * sets,steps-Ki_start+1);%%same with orstd
%     save(sprintf('highres_kura_test_Fg=%.3f.mat', freq_gap(FG)), 'Or_std_data_reshaped');
end




figure();
for fig = 1:size(Or_data_reshaped,1)
    plot(K_vec(Ki_start:end), Or_data_reshaped(fig,:))
    hold on
end

