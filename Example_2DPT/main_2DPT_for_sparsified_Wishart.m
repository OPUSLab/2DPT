clear all
close all
clc


%% 2D-PT for sparsified Wishart instances
% Corentin Delacour, OPUSlab, University of California, Santa Barbara
% delacour@ucsb.edu

% Add path for Wishart dataset
addpath('../Wishart_instances/')

% Add path for source code
addpath('../Source_code')

% Number of logical nodes in original graph (40-100)
N=30

% Number of copy nodes per original node
N_copy=2

% Total number of physical nodes
Nt=N_copy*N

% Instance index
instance=1

% Load Wishart instance
alpha=0.75; % Wishart hardness parameter
file_J='../Wishart_instances/instances_alpha'+string(alpha)+'/J_N'+string(N)+'.mat';
file_g='../Wishart_instances/instances_alpha'+string(alpha)+'/ground_N'+string(N)+'.mat';
[J,ground_energy]=load_Wishart_instance(N,instance,file_J,file_g);
h=zeros(N,1); % No external field


%% Running Adaptive algorithm (APT)
Mb=20; % max number of beta values
Nw=20; % max number of W0 values

% Initial values for APT
W0_init=N_copy*N/5;
beta_init=0.1;

Nchain=20; % Number of Markov chains for Variance estimation

% Number of MCS for variance estimation
if N<=100
    MCS=4e3; 
else
    MCS=8e3;
end

% "Learning rate" or "step" for APT
alpha_beta=1;
alpha_W0=0.4;

% Start parallel pool for APT
p = gcp("nocreate");
if isempty(p)
    parpool(10)
end

%Running APT
[beta,W0] = APT_2D_for_sparsified_graph(beta_init,W0_init,alpha_beta,alpha_W0,Nchain,MCS,Mb,Nw,J,h,N_copy,Nt);

Nw=length(W0)
Mb=length(beta)

figure
subplot(1,2,1)
plot(1:Nw,W0,'--o','LineWidth',2)
xlabel('Column index')
ylabel('Copy strength W0')
subplot(1,2,2)
plot(1:Mb,beta,'--o','LineWidth',2)
xlabel('Row index')
ylabel('\beta')
fontsize(20,"points")



%% Running 2D PT with final profile
% Number of Monte Carlo sweeps
MCS=1e4;

% Sweep-to-swap ratio
MCS_per_swap=5;

% Number of swaps
Nswap=floor(MCS/MCS_per_swap); 


Pswap_b=zeros([Nw,Mb,Mb]);
Pswap_p=zeros([Mb,Nw,Nw]);

S_init=sign(randn([length(beta) length(W0) Nt])); % Markov chain start

% Running 2D-PT
[spin_store_swap,count_swap_b,count_swap_p,energy,constraint_g,copy_idx]  = dynamics_2D_PT_sparsified_graph(N,N_copy,J,h,beta,W0,MCS_per_swap,MCS,S_init);


% Reading min energy from bottom-right replica that satisfies constraints
% (=feasible)
min_meas_energy=min(energy(:,end,end).*double(constraint_g(:,end,end)==0),[],'all');

% Computing back the original Wishart energy
min_original_energy=min_meas_energy+N*(N_copy-1)*W0(end)

% Wishart ground energy
ground_energy

% Normalized residual energy
norm_residual_energy=(min_original_energy-ground_energy)/abs(ground_energy)


% Plotting energy for each column
figure
for column=1:Nw
    subplot(1,Nw,column)
    leg=[];
    for t=1:length(beta)
        last_column_energy=squeeze(energy(:,:,column));
        plot(1:Nswap,last_column_energy(:,t));
        hold on
        leg=[leg "\beta="+string(beta(t))];
    end
    legend(leg)
    xlabel('Swap')
    ylabel('Energy')
    %ylim([min_energy max_energy ])
    title('Column '+string(column))
end


% Plotting constraints for each row
figure
for row=1:Mb
    leg=[];
    subplot(1,Mb,row)
    for p=1:length(W0)
        last_column_g=squeeze(constraint_g(:,row,:));
        plot(1:Nswap,last_column_g(:,p));
        hold on
        leg=[leg "W0="+string(W0(p))];
    end
    legend(leg)
    xlabel('Swap')
    ylabel('Constraint g')
    title('Row '+string(row))
end


% Swap probabilities along penalties
for t=1:length(beta)
    for p=1:Nw-1
        Pswap_p(t,p,p+1)=4*count_swap_p(t,p,p+1)/Nswap;
    end
end

% Swap probabilities along beta values
for p=1:length(W0)
    for t=1:Mb-1
        Pswap_b(p,t,t+1)=4*count_swap_b(p,t,t+1)/Nswap;
    end
end

%% Plot swap probabilities

% For plot purposes
cols_beta = ceil(sqrt(Mb));
rows_beta = ceil(Mb / cols_beta);

cols_W0 = ceil(sqrt(Nw));
rows_W0= ceil(Nw / cols_W0);


figure
for t=1:length(beta)
subplot(rows_beta,cols_beta,t)
swap_matrix=[];
lab=[];
for p=1:Nw-1
    swap_matrix=[swap_matrix ; Pswap_p(t,p,p+1)];
    lab=[lab string(p)+'-'+string(p+1)];
end
bar([1:1:Nw-1],swap_matrix)
xticklabels(lab)
ylabel('Swap probability')
xlabel('Replica pair')
title('row '+string(t))
ylim([0 1])
grid on
end

figure
for p=1:length(W0)
subplot(rows_W0,cols_W0,p)
swap_matrix=[];
lab=[];
for t=1:Mb-1
    swap_matrix=[swap_matrix ; Pswap_b(p,t,t+1)];
    lab=[lab string(t)+'-'+string(t+1)];
end
bar([1:1:Mb-1],swap_matrix)
grid on
xticklabels(lab)
ylabel('Swap probability')
xlabel('Replica pair')
title('column '+string(p))
ylim([0 1])
end



