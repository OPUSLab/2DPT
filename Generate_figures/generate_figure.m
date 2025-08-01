%% Generate figures from Main manuscript: "Two-dimensional Parallel Tempering for Constrained Optimization"
% Corentin Delacour et al., OPUS lab, University of California, Santa
% Barbara


clear all
close all
clc


%%%%%%%%%%%%%%%%%%% FIG.2 %%%%%%%%%%%%%%%%%%%%%%%%%%

data_R4_noswap=readmatrix('../Processed_results/data_FA_noswap_R4.txt');
data_R4_swap=readmatrix('../Processed_results/data_FA_swap_R4.txt');
    
MCS_array=data_R4_swap(:,1);
KL_original_noswap_m4=data_R4_noswap(:,2);
ci_original_no4=data_R4_noswap(:,3:4);

KL_original_swap_m4=data_R4_swap(:,2);
ci_original_4=data_R4_swap(:,3:4);

figure(1)
shadedErrorBar(MCS_array,KL_original_noswap_m4,([(ci_original_no4(:,2)-KL_original_noswap_m4) (KL_original_noswap_m4-ci_original_no4(:,1))])','lineProps',{'-o','color',[0 0.4470 0.7410]})
grid on
hold on
shadedErrorBar(MCS_array,KL_original_swap_m4,([(ci_original_4(:,2)-KL_original_swap_m4) (KL_original_swap_m4-ci_original_4(:,1))])','lineProps',{'-o','color',[0.8500 0.3250 0.0980]})
xlabel('Monte Carlo Sweeps')
ylabel('KL Divergence')
xscale('log')
yscale('log')
legend('No swaps','With swaps')
title('Sampling a sparsified Full Adder')
    




%%%%%%%%%%%%%%%%%%%%% FIG. 3 %%%%%%%%%%%%%%%%%%%%%%%%% 

%% Plot minimum energies for 2D-PT and J-column PT

colors_plot_1D=[198,219,239
158,202,225
107,174,214
66,146,198
33,113,181
8,81,156
8,48,107]/255;

colors_plot_2D=[252,187,161
252,146,114
251,106,74
239,59,44
203,24,29
165,15,21
103,0,13]/255;

index_size=0;
for N=40:10:100

    index_size=index_size+1;
    data_1D=readmatrix('../Processed_results/data_energy_1D_N'+string(N)+'.txt');
    MCS=data_1D(:,1);
    energy_1D=data_1D(:,2);
    energy_ci_low_1D=data_1D(:,3);
    energy_ci_high_1D=data_1D(:,4);

    figure(2)
    subplot(1,2,2)
    shadedErrorBar(MCS,energy_1D,[(energy_ci_high_1D-energy_1D) (energy_1D-energy_ci_low_1D) ],'lineProps',{'-o','color',colors_plot_1D(index_size,:)})
    hold on
    grid on
    xscale('log')
    yscale('log')
    xlabel('Monte Carlo Sweeps')
    ylabel('Residual Energy \rho_E')
    title('J-column PT')

    data_2D=readmatrix('../Processed_results/data_energy_2D_N'+string(N)+'.txt');
    MCS=data_2D(:,1);
    energy_2D=data_2D(:,2);
    energy_ci_low_2D=data_2D(:,3);
    energy_ci_high_2D=data_2D(:,4);

    subplot(1,2,1)
    shadedErrorBar(MCS,energy_2D,[(energy_ci_high_2D-energy_2D) (energy_2D-energy_ci_low_2D) ],'lineProps',{'-o','color',colors_plot_2D(index_size,:)})
    hold on
    grid on
    xscale('log')
    yscale('log')
    xlabel('Monte Carlo Sweeps')
    ylabel('Residual Energy \rho_E')
    title('2D-PT')

end
subplot(1,2,1)
legend('N=40','N=50','N=60','N=70','N=80','N=90','N=100')
subplot(1,2,2)
legend('N=40','N=50','N=60','N=70','N=80','N=90','N=100')

%% Collapsing data and plotting the results

% Starting data point for collapse
start_point=7

% J-column PT Collapse

% Dynamic exponent
mu=11

% Static exponent
b=0

figure(3)
subplot(1,2,2)
size_index=0;
for N=50:10:100

size_index=size_index+1;

data_1D=load('../Processed_results/data_energy_1D_N'+string(N)+'.txt');
MCS=data_1D(start_point:end,1);
energy_1D=data_1D(start_point:end,2);

loglog(MCS*N^(-mu),energy_1D*N^b,'o','LineWidth',2,'Color',colors_plot_1D(size_index,:),'MarkerFaceColor',colors_plot_1D(size_index,:))
hold on

end
xlabel('t N^{-\mu}')
ylabel('Residual Energy')
title('N-column PT: \mu='+string(mu))
legend('N=50','N=60','N=70','N=80','N=90','N=100')


% 2D-PT Collapse

% Dynamic exponent
mu=6

% Static exponent
b=0

subplot(1,2,1)


size_index=0;
for N=50:10:100

size_index=size_index+1;
data_2D=load('../Processed_results/data_energy_2D_N'+string(N)+'.txt');
MCS=data_2D(start_point:end,1);
energy_2D=data_2D(start_point:end,2);

loglog(MCS*N^(-mu),energy_2D*N^b,'o','LineWidth',2,'Color',colors_plot_2D(size_index,:),'MarkerFaceColor',colors_plot_2D(size_index,:))
hold on

end
xlabel('t N^{-\mu}')
ylabel('Residual Energy')
title('2D-PT: \mu='+string(mu))
legend('N=50','N=60','N=70','N=80','N=90','N=100')
