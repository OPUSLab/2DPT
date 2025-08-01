function [J,ground_energy] = load_Wishart_instance(N,instance_number,file_J,file_g)
% load Wishart instance with parameter alpha


Jdat=load(file_J,"J");

Jall=Jdat.J;
J=-N*squeeze(Jall(instance_number,:,:)); % factor N to counter the normalization done during generation, minus sign to counter positive Hamiltonian

ground_energy_dat=load(file_g,"ground");

ground_energy=N*ground_energy_dat.ground(instance_number); % factor N to counter the normalization done during generation

end
