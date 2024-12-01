%%%CODE SAMPLE. THIS SCRIPT WAS USED DURING ATTEMPTS TO UTILIZE UofM's SUPERCOMPUTING CLUSTER FOR COMPUTATIONALLY INTENSIVE KURAMOTO MODEL SIMULATIONS. 
%%%THE SCRIPT ENABLES OPTIMAL USAGE OF AN INFINITE NUMBER OF CORES FOR AN INFINITE NUMBER OF COUPLING STRENGTHS IN DISCONTINUOUS KURAMOTO MODEL SIMULATIONS


function kura_parallel_supercomputer(Kura_gen_values, savepath, sigSavepath)
tic;
w=load(sprintf('NaturalFrequencies_FreqGap=%.3f', Kura_gen_values(2))); %%load all iterations of natural frequencies
conn_matrices = load(sprintf('MAT_targetFG=%.3f',Kura_gen_values(2)));
MAT = conn_matrices.MAT(:,:,Kura_gen_values(3));
MAT = Kura_gen_values(1)*MAT;
N = size(MAT,1);

T=0:0.01:80;
dt=T(2)-T(1);
transient_time=40;
noise_std=0;
sig=zeros(N,length(T)-transient_time*1/dt);
noisein=noise_std*randn(N,length(T));


start = 2*pi*rand(N,1)-pi;


 [phi] = Kuramoto_noisein(MAT,T,w.W_save(Kura_gen_values(3),:)',noisein,start);
    phi=phi(:,transient_time*1/dt+1:end);
    [or,op,or_std,or_t] = OrderParameter2(sin(phi)');
    Or=or;
    Or_std=or_std;

    sig(:,:)=phi;
    
    t=T(1:5:end-transient_time*1/dt);
    sig=sin(sig(:,:)); 
    fname=sprintf('%sKura_discontinuous_supercomputer_k=%.6f_freq_gap=%.3f_iter=%.i.mat', ...
        savepath, Kura_gen_values(1), Kura_gen_values(2), Kura_gen_values(3) );
    fname2=sprintf('%sKura_discontinuous_supercomputer_sig_k=%.6f_freq_gap=%.3f_iter=%.i.mat',sigSavepath, ...
        Kura_gen_values(1),Kura_gen_values(2),Kura_gen_values(3));
    save(fname,'t','Or','Or_std')
    save(fname2,'sig')
    toc;
