%%CODE SAMPLE. THIS GENERATES TIME SERIES AND VARIANCE IN ORDER PARAMETER FOR LARGE NUMBER OF KURAMOTO MODEL SIMULATIONS.

%%first ripser test 4.5

clear; clc;

addpath('D:\Joe\Zauberbaum\AdjustFrequencyGapNetworks\ripser_test_4.5\')
target_FG = .35;
K_vec = [linspace(0,.15,2000)]; %%%k resolution


%%savepaths
savepath= 'D:\Joe\Zauberbaum\first_ripser_test\ripser_test_4.5_or_orstd\';
Sig_savepath= 'D:\Joe\Zauberbaum\first_ripser_test\ripser_test_4.5_signal\';
filedir = savepath;
addpath(Sig_savepath);

for i = 1:length(target_FG)
load(sprintf('NaturalFrequencies_and_MAT_FreqGap=%.3f.mat', target_FG(i)));
disp(['running simulation for target FG ' num2str(target_FG(i))])

sets = 5; 

T=0:0.01:80;
dt=T(2)-T(1);
transient_time=40; cores = 16; %% 16 cores
iterations = sets * cores;
tic
noise_std=0; %%%no noise in 2013 paper we are following 10.12

or =zeros(length(K_vec),iterations);
or_std = zeros(length(K_vec),iterations);
K_vec_length = length(K_vec);



for set = 1:sets
parfor iter = 1 + ((set-1)*cores): set * cores
    w = W_save(iter,:);
    MAT0 = MAT(:,:,iter);
    N=size(MAT,1);
    sig=zeros(N,length(T)-transient_time*1/dt);
% for Ki=1:K_vec_length
 for Ki = 1:K_vec_length %%%first ripser test 3.29
    k=K_vec(Ki);
    MAT2=k*MAT0;
    % noise
    noisein=noise_std*randn(N,length(T));
    if Ki ==1
   start = 2*pi*rand(N,1)-pi;
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
        savepath, target_FG(i), k, iter, Ki);
     fname2=sprintf('%sKura_Sig_freq_gap=%.4f_k=%.6f_iteration=%.3fstate=%.2f.mat',Sig_savepath,target_FG(i),k,iter,Ki)
     parsave_kura_data(fname,t,w,k,MAT2,Or_std_vec,Or_vec)
     parsave_kura_data_signal(fname2,sig)
    
end
end
toc
end
end
