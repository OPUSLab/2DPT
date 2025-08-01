function [spin_store,current_e,current_g] = dynamics_single_replica_sparsification3(N,Nt,Js,hs,copy_idx,beta,MCS,S_init)
% N: number of spins in initial graph
% Nt: total number of spins after sparsification
% Js: J matrices for each all replicas (1st and 2nd dimensions: beta and
% penalty, respectively)
% hs: h vectors for all replicas (1st and 2nd dimensions: beta and
% penalty, respectively)
% beta: vector of beta values
% W0: vector of copy values
% copy_idx: array of copy indices (row i is original spin)
% Nswap: number of exchanges in both beta and P directions
% MCS: total number of Monte Carlo sweeps

kmax=MCS+1;

% including init
current_g=zeros([MCS+1 1]);
current_e=zeros([MCS+1 1]);

spin_store=zeros(MCS+1, Nt );



current_J=reshape(Js,[Nt Nt]);
%current_J2=reshape(Js2(t,p,:,:),[Nt Nt]);
current_h=reshape(hs,[Nt 1]);
S_temp_pt=S_init;
I1=zeros(Nt,1);
rand_matrix=rand([Nt,kmax]);

for k=1:kmax
    spin_store(k,:)=S_temp_pt;

    %Chain breaking g calculation
    current_copy_idx=copy_idx;
    sz_copy_idx=size(current_copy_idx);
    Ncopy=sz_copy_idx(2);
    if Ncopy>=N
        fprintf('Problem with copy index matrix\n')
    end
    for i=1:N
        j=1;
        current_copy=current_copy_idx(i,1);
        previous_copy=i;
        while (j<=Ncopy)&&(current_copy~=0)
            current_g(k)=current_g(k)+0.25*(S_temp_pt(previous_copy)-S_temp_pt(current_copy))^2; % penalty function for sparsification
            j=j+1;
            previous_copy=current_copy;
            if j<=Ncopy
                current_copy=current_copy_idx(i,j);
            end
        end
    end
    % Energy calculation
    current_e(k)=-0.5*S_temp_pt'*Js*S_temp_pt-hs*S_temp_pt;%+extra;%+constant;

        %p-bit sweep
    for i=1:Nt
        I1(i)=beta*(current_J(i,:)*S_temp_pt+current_h(i)); % input current without copy connections
        %I2(i)=current_J2(i,:)*S_temp_pt; % only copies
        S_temp_pt(i)=sign(tanh(I1(i))-2*rand_matrix(i,k)+1); % spin update
    end
    %state_store(n,t,p,k)=0.5*(reshape(S_temp(t,p,:),[1 Nt])+1)*(binary_weights)+1;%binaryVectorToDecimal((reshape(S_temp(t,p,:),[1 Nt])+1)/2)+1;
   
             
end

        

end