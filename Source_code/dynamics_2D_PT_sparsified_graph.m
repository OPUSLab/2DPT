function [spin_store_swap,count_swap_b,count_swap_p,current_e,current_g,copy_idx] = dynamics_2D_PT_sparsified_graph(N,N_copy,J,h,beta,W0,MCS_per_swap,MCS,S_init)
% 2D-PT algorithm for sparsified graphs
% Corentin Delacour, OPUSlab, University of California, Santa Barbara
% delacour@ucsb.edu

% INPUTS:
% N: number of spins in initial graph
% N_copy: number of copy spins per original spin
% J and h: Ising matrix and external field
% beta: vector of beta values
% W0: vector of copy values
% MCS: total number of Monte Carlo sweeps
% MCS_per_swap: number of MCS between two swaps
% S_init: initial spin values of dimension [Mb Nw Nt]

% OUTPUTS:
% spin_store_swap: spin states after each swap. Dim [Nswap Mb Nw Nt]
% count_swap_b/p: count swaps along beta/penalty axis for every replica pair
% Dim [Nw Mb Mb] and [Mb Nw Nw]
% current_e: energy at every swap
% current_g: constraint function g at every swap (measures copy
% disagreement)


% Sparsifying graph for each penalty strength (copy value W0)
for p=1:length(W0)
    [Js(p,:,:), hs(p,:),copy_idx(p,:,:),~]=sparse_J_h_input_copies(J,h, W0(p), N_copy);
end

% total number of spins
Nt=length(hs(1,:));


count_swap_p=zeros([length(beta),length(W0),length(W0)]);
count_swap_b=zeros([length(W0),length(beta),length(beta)]);


%Max number of Monte Carlo sweeps before swap
kmax=MCS_per_swap;

Nswap=floor(MCS/MCS_per_swap); % Total number of swaps

S_temp=S_init;

spin_store_swap=zeros([Nswap length(beta) length(W0) Nt]);
current_g=zeros([Nswap length(beta) length(W0)]);
current_e=zeros([Nswap length(beta) length(W0)]);

sweep_index=zeros(length(beta), length(W0));

% Local inputs for spins
I1=zeros(Nt,1);

% to alternate swap direction
swap_direction=-1;

for n=1:Nswap
    for t=1:length(beta)
        for p=1:length(W0)
            rand_matrix=rand(Nt,kmax);
            current_J=reshape(Js(p,:,:),[Nt Nt]);
            %current_J2=reshape(Js2(t,p,:,:),[Nt Nt]);
            current_h=reshape(hs(p,:),[Nt 1]);
            S_temp_pt=reshape(S_temp(t,p,:),[Nt 1]);
            for k=1:kmax
                sweep_index(t,p)=sweep_index(t,p)+1;
                %p-bit sweep
                for i=1:Nt
                    I1(i)=beta(t)*(current_J(i,:)*S_temp_pt+current_h(i)); % input current without copy connections
                    %I2(i)=current_J2(i,:)*S_temp_pt; % only copies
                    S_temp_pt(i)=sign(tanh(I1(i))-2*rand_matrix(i,k)+1); % spin update
                end
                %state_store(n,t,p,k)=0.5*(reshape(S_temp(t,p,:),[1 Nt])+1)*(binary_weights)+1;%binaryVectorToDecimal((reshape(S_temp(t,p,:),[1 Nt])+1)/2)+1;
                %spin_store(t,p,sweep_index(t,p),:)=S_temp_pt;
            end
            S_temp(t,p,:)=S_temp_pt;
           
        end
    end
    

    % energy for each replica before swap
    for t=1:length(beta)
            for p=1:length(W0)
                %constant=W0(p)*2*(Nt-N);% P(S1-S2)^2 creates constant P(S1^2+S2^2)
                S_temp_pt=reshape(S_temp(t,p,:),[Nt 1]);
                current_copy_idx=reshape(copy_idx(p,:,:),N,[]);
                sz_copy_idx=size(current_copy_idx);
                Nc=sz_copy_idx(2);
                if Nc>=N
                    fprintf('Problem with copy index matrix\n')
                end
                for i=1:N
                    j=1;
                    current_copy=current_copy_idx(i,1);
                    previous_copy=i;
                    while (j<=Nc)&&(current_copy~=0)
                        
                        current_g(n,t,p)=current_g(n,t,p)+0.25*(S_temp_pt(previous_copy)-S_temp_pt(current_copy))^2; % penalty function for sparsification
                        j=j+1;
                        previous_copy=current_copy;
                        if j<=Nc
                            current_copy=current_copy_idx(i,j);
                        end
                    end
                end
    
                %extra=W0(p)*(N-current_g(n,t,p))-current_g(n,t,p)*W0(p);
           
                current_e(n,t,p)=-0.5*S_temp_pt'*reshape(Js(p,:,:),[Nt Nt])*S_temp_pt-S_temp_pt'*reshape(hs(p,:),[Nt 1]);%+extra;%+constant;
             
            end
    end

    if mod(n,2)==0
        even_swap=1;
    else
        even_swap=0;
    end

    if swap_direction==1
        %swap probabilities in Penalty direction
        prob_row_swap=zeros(length(beta),length(W0),length(W0));
        for t=1:1:length(beta)
            % Penalty swap (fixed row)
            for l=1+even_swap:2:length(W0)-1
                dW=(W0(l+1)-W0(l));
                prob_row_swap(t,l,l+1)=min([1 exp(beta(t)*(current_g(n,t,l+1)-current_g(n,t,l))*dW)]);% Metropolis rule
                % Glauber also possible: 1/(1+exp(-beta(t)*(current_g(n,t,l+1)-current_g(n,t,l))*dW)); 
                if rand(1)<= prob_row_swap(t,l,l+1) %swap
                    count_swap_p(t,l,l+1)=count_swap_p(t,l,l+1)+1;
                    S_buffer=S_temp(t,l+1,:);
                    S_temp(t,l+1,:)=S_temp(t,l,:);
                    S_temp(t,l,:)=S_buffer;
                end
             end
        end
    end

    if swap_direction==-1
        %swap probabilities in temperature direction
        prob_column_swap=zeros(length(W0),length(beta),length(beta));
        for p=1:length(W0)
            % beta swap (fixed column)
            for l=1+even_swap:2:length(beta)-1
                dbeta=beta(l+1)-beta(l); % beta_large-beta_small
                prob_column_swap(p,l,l+1)=min([1 exp(dbeta*(current_e(n,l+1,p)-current_e(n,l,p)))]);% Metropolis rule
                % Glauber also possible: 1/(1+exp(-dbeta*(current_e(n,l+1,p)-current_e(n,l,p))));
                if rand(1)<= prob_column_swap(p,l,l+1) %swap
                    count_swap_b(p,l,l+1)=count_swap_b(p,l,l+1)+1;
                    S_buffer=S_temp(l,p,:);
                    S_temp(l,p,:)=S_temp(l+1,p,:);
                    S_temp(l+1,p,:)=S_buffer;
                    %fprintf("beta swap\n")
                end
                   
            end
        end
    end


    for t=1:length(beta)
        for p=1:length(W0)
            % Storing spins at each swap
            spin_store_swap(n,t,p,:)=squeeze(S_temp(t,p,:));
        end
    end

    % Alternate swap direction
    if mod(n,2)==0
        swap_direction=-swap_direction;
    end

end

end