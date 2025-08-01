function [beta,W0] = APT_2D_for_sparsified_graph(beta_init,W0_init,alpha_beta,alpha_W0,Nchain,MCS,Mb,Nw,J,h,N_copy,Nt)
% Adaptive schedule for 2D-PT for sparsified graphs with copy edge W0=penalty strength
% Corentin Delacour, OPUSlab, University of California, Santa Barbara
% delacour@ucsb.edu

% This function outputs the beta and W0 values that are computed
% iteratively from initial beta_init and W0_init values
% J and h are the Ising matrix and external field
% N_copy is the number of copies per original node
% Nt is the total number of spins after sparsification
% Nchain is the number of chain per replica used to estimate the standard
% deviation of energy and constraint g
% Mb is the max number of beta values
% Nw is the max number of W0 values
% alpha parameters are the learning rates or increments

% Initial arrays that will be reshaped at the end
W0=ones(1,Nw)*W0_init;
beta=ones(1,Mb)*beta_init;

N=length(h);


%% getting the number of beta values for the first column

% Graph sparsification
[Js, hs, cpy_idx, ~] = sparse_J_h_input_copies(J,h, W0(1), N_copy);

% Min standard deviation serving as threshold to stop adaptive schedule
min_std=mean(abs(Js(Js~=0)),'all')/50; %/10

% stored energy at swaps
E=zeros(Mb,Nchain);
g=zeros(Nw,Nchain);

S_final=sign(randn([round(Nt),Nchain])); % Markov chain start for the column

n=1;
m=1;
while m<=Mb % for each beta value (row)
    dE=0;
    dg=0;
    mg=0;
        parfor c=1:Nchain
            S_init=S_final(:,c)
            % Running replica
            [spin_store,current_e,current_g] = dynamics_single_replica_sparsified(N,Nt,Js,hs,cpy_idx,beta(m),MCS,S_init);
            % continuing Markov chain for adjacent beta-replicas and
            % adjacent columns
            S_final(:,c)=reshape(spin_store(end,:),[Nt,1]);

            if (m==1) % removing chain initialization or "burn-in" period
                current_e=current_e(round(MCS/2):end);
                current_g=current_g(round(MCS/2):end);
            end
            dE=dE+std(current_e);
            dg=dg+std(current_g);
            mg=mg+mean(current_g);
          
        end

        % population average
        dE=dE/Nchain;
        dg=dg/Nchain;
        mg=mg/Nchain;

       % figure
       %  subplot(2,1,1)
       %  plot(current_g,'LineWidth',2)
       %  xlabel('Sample')
       %  ylabel('g')
       %  title('\sigma_g='+string(dg(m)))
       %  fontsize(20,"points")
       %  ylim([0 40])
       %  subplot(2,1,2)
       %  plot(current_e,'LineWidth',2)
       %  xlabel('Sample')
       %  ylabel('E')
       %  title('\sigma_E='+string(dE(m)))
       %  fontsize(20,"points")

        if (m<=Mb-1) && (dE>min_std)
            W0_store(m,2)=W0(1)+alpha_W0/(beta(m)*dg); % storing W0 values computed for this beta value
            beta(m+1)=beta(m)+alpha_beta/dE;
            m=m+1;
        else
            break;
        end


end

% final number of rows
Mb=m;
beta=beta(1:Mb);
beta_store(:,1)=beta;


W0(2)=median(W0_store(:,2));
n=2;
while n<=Nw % for each W0 value (column)
 
    S_final=sign(randn([round(Nt),Nchain])); % Markov chain start for the column
    
    % Sparsifying with current W0 value
    [Js, hs, cpy_idx, ~] = sparse_J_h_input_copies(J,h, W0(n), N_copy);

    dE=zeros(Mb,1);
    dg=zeros(Mb,1);
    mg=zeros(Mb,1);


    for m=1:Mb % for each beta value (row)
        dE=0;
        dg=0;
        mg=0;
        parfor c=1:Nchain
            S_init=S_final(:,c)

            % Running replica
            [spin_store,current_e,current_g] = dynamics_single_replica_sparsified(N,Nt,Js,hs,cpy_idx,beta(m),MCS,S_init);

            % continuing Markov chain for adjacent beta-replicas and
            % adjacent columns
            S_final(:,c)=reshape(spin_store(end,:),[Nt,1]);

            if (n==1)&(m==1) % removing chain initialization
                current_e=current_e(round(MCS/2):end);
                current_g=current_g(round(MCS/2):end);
            end
            dE=dE+std(current_e);
            dg=dg+std(current_g);
            mg=mg+mean(current_g);

        end
        % population average
        dE=dE/Nchain;
        dg=dg/Nchain;
        mg=mg/Nchain;

        if m<=Mb-1
            if dE>min_std
                beta(m+1)=beta(m)+alpha_beta/dE;
            else
                % Thresold is reached but we cannot stop here: we need to
                % complete the rows
                beta(m+1)=2*beta(m); % arbitrary factor: geometric suite
            end

        end
        if n<=Nw-1
            if dg>min_std
                W0_store(m,n+1)=W0(n)+alpha_W0/(beta(m)*dg); % storing W0 values computed for this beta value
            else
                % Thresold is reached but we cannot stop here: we need to
                % complete the rows
                W0_store(m,n+1)=2*W0(n); % arbitrary factor: geometric suite
            end
        end

    end

    beta_store(:,n)=beta;


    if mg(end)<0.000001 % checking when mean constraint is near 0
        fprintf('last replica from last column achieves near g=0\n')
        break
    else
        if n<=Nw-1
             W0(n+1)=median(W0_store(:,n+1));
        end
        n=n+1;
    end


end

 % taking median value of beta values
 beta=median(beta_store,2)';

 if n<Nw
    W0=W0(1:n);
 end

end
 