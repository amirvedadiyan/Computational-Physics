clear; 
clc; 
close all;

%-------initialization-------
Nx = 4;
Ny = 4;
J1 = 1; % Nearest neighbor Heisenberg coupling
J2 = 0.3; % Next nearest neighbor Heisenberg coupling
bnd = 1; % boundary condition : 1--> PBC  0 --> OBC  -1 --> APBC
hbar = 1;

%-------pauli matrices-------
sigma_0 = speye(2);
sigma_z = sparse([1,0; 0,-1]);
sigma_x = sparse([0, 1; 1,0]);
sigma_y = sparse([0,-1i; 1i,0]);
sigma_p = (sigma_x + 1i * sigma_y)/2;
sigma_m = (sigma_x - 1i * sigma_y)/2;
n = (sigma_p * sigma_m); % number operator

%-------general problem-> spin operators for n sites-------
Sx{1} = hbar/2 * sigma_x;
Sy{1} = hbar/2 * sigma_y;
Sz{1} = hbar/2 * sigma_z;
ni{1} = hbar/2 * n;
id = sigma_0;
for num_sites = 2:Nx*Ny
    for pos = 1:num_sites-1
        Sx{pos} = kron(Sx{pos}, sigma_0);
        Sy{pos} = kron(Sy{pos}, sigma_0);
        Sz{pos} = kron(Sz{pos}, sigma_0);
        ni{pos} = kron(ni{pos}, sigma_0);
    end
    Sx{num_sites} = kron(id, hbar/2*sigma_x);
    Sy{num_sites} = kron(id, hbar/2*sigma_y);
    Sz{num_sites} = kron(id, hbar/2*sigma_z);
    ni{num_sites} = kron(id,hbar/2*n);
    id = kron(id, sigma_0);
end
%-------Hamiltonian-------
H = sparse(0);
for column = 1:Nx
    for row = 1:Ny
        index = (column-1)*Ny + row;
        % nearest neighbors interaction
        if row == Ny && column == Nx 
            % periodic bc
            H = H + J1 * (Sx{index}*Sx{index-Ny+1}+Sy{index}*Sy{index-Ny+1}+Sy{index}*Sy{index-Ny+1});
            H = H + J1 * (Sx{index}*Sx{Ny}+Sy{index}*Sy{Ny}+Sy{index}*Sy{Ny});
        elseif row == Ny
            % only right + periodic bc
            H = H + J1 * (Sx{index}*Sx{index+Ny}+Sy{index}*Sy{index+Ny}+Sz{index}*Sz{index+Ny});
            H = H + J1 * (Sx{index}*Sx{index-Ny+1}+Sy{index}*Sy{index-Ny+1}+Sz{index}*Sz{index-Ny+1});
        elseif column == Nx
            % only up + periodic bc
            H = H + J1 *(Sx{index}*Sx{index+1}+Sy{index}*Sy{index+1}+Sz{index}*Sz{index+1});
            H = H + J1 * (Sx{index}*Sx{row}+Sy{index}*Sy{row}+Sz{index}*Sz{row});
        else
            % right and up
            H = H + J1 * (Sx{index}*Sx{index+Ny}+Sy{index}*Sy{index+Ny}+Sz{index}*Sz{index+Ny}) + J1 *(Sx{index}*Sx{index+1}+Sy{index}*Sy{index+1}+Sz{index}*Sz{index+1});
        end
        % next nearest neighnors interaction
        if column == Nx
            % periodic bc for J2 
            if row == 1 
                H = H + J2 * (Sx{index}*Sx{row+1}+Sy{index}*Sy{row+1}+Sz{index}*Sz{row+1});
            elseif row == Ny
                H = H + J2 * (Sx{index}*Sx{row-1}+Sy{index}*Sy{row-1}+Sz{index}*Sz{row-1});
            else 
                H = H + J2 * (Sx{index}*Sx{row+1}+Sy{index}*Sy{row+1}+Sz{index}*Sz{row+1}) + (Sx{index}*Sx{row-1}+Sy{index}*Sy{row-1}+Sz{index}*Sz{row-1});
            end
        elseif row == Ny 
            % only down right
            H = H + J2 * (Sx{index}*Sx{index+Ny-1}+Sy{index}*Sy{index+Ny-1}+Sz{index}*Sz{index+Ny-1});
            % periodic bc for J2
            if column == 1
                H = H + J2 * (Sx{index}*Sx{index+1}+Sy{index}*Sy{index+1}+Sz{index}*Sz{index+1});
            elseif column == Nx
                H = H + J2 * (Sx{index}*Sx{index-2*Ny+1}+Sy{index}*Sy{index-2*Ny+1}+Sz{index}*Sz{index-2*Ny+1});
            else
                H = H + J2 * (Sx{index}*Sx{index+1}+Sy{index}*Sy{index+1}+Sz{index}*Sz{index+1}) + J2 * (Sx{index}*Sx{index-2*Ny+1}+Sy{index}*Sy{index-2*Ny+1}+Sz{index}*Sz{index-2*Ny+1});
            end
        elseif row == 1
            % only up right 
            H = H + J2 * (Sx{index}*Sx{index+Ny+1}+Sy{index}*Sy{index+Ny+1}+Sz{index}*Sz{index+Ny+1});
        else 
            % both up right and down right
            H = H + J2 * (Sx{index}*Sx{index+Ny-1}+Sy{index}*Sy{index+Ny-1}+Sz{index}*Sz{index+Ny-1}) + J2 * (Sx{index}*Sx{index+Ny+1}+Sy{index}*Sy{index+Ny+1}+Sz{index}*Sz{index+Ny+1});
        end
   end
end

%-------total number operator-------
N_tot = sparse(0);
L = Nx * Ny;
for i =1:L
    N_tot = N_tot + ni{i};
end
N_tot_diag = diag(N_tot);
%States with Sz{tot} = 0 is N_tot = L/2*h/2
F = N_tot_diag==(L*hbar/4);
H_tilde = H(F,F);
[eigenvectors,eigenvalues] = eigs(H_tilde,size(H_tilde,1),'sa');
E = diag(eigenvalues);

%plot

figure(1);
nE = min(20,size(H,1));
a = zeros(nE,1);
b = E(1:nE);
plot(a,b, '--o','MarkerSize',10,...
    'MarkerEdgeColor','black',...
    'MarkerFaceColor', [0.5 0.5 0.5])
set(gca,'FontSize',16)
ylabel('E(n)');
grid on;
box on;
title('Lowest 20 energy eigenvalues');




u = sparse([], [], [], 2^L, 2^L);
u(F,F) = eigenvectors;
[row_indices, ~] = find(F);
% Extract the first row index
first_psi_index = min(row_indices);
psi = u(:,first_psi_index);
%--------correlation function-------- ***change this for x,y direction or
%number difference***
SS_corr = zeros(Ny);
mc = zeros(Ny);
nn = 1;
for i1=1:Ny:L
     O = Sx{1}*Sx{i1}+Sy{1}*Sy{i1}+Sz{1}*Sz{i1};
     mc(nn) = psi' * O * psi; 
     nn = nn +1;
end
for i1 = 1:Ny:L
    for i2 = i1:Ny:L 
        O = Sx{i1}*Sx{i2}+Sy{i1}*Sy{i2}+Sz{i1}*Sz{i2};
        SS_corr(i1,i2) = psi' * O * psi; 
        SS_corr(i2,i1) = SS_corr(i1,i2);
    end
end
figure(4);
x = (0:Ny-1);
y = mc(:,1);
plot(x,y,'-d');
xlabel('distance in x direction')
ylabel('Spin-Spin correlation')



L = Nx*Ny;
ens = zeros(L,1);
for dim_L = 0:L
    dim_R = L-dim_L;
    psi_tilde = reshape(psi,2^dim_R,2^dim_L);
    rho_L = conj(psi_tilde'*psi_tilde);
    [V,E] = eigs(rho_L,size(rho_L,1),'sa');
    E(E < 1e-13) = 1e-13;
    E = diag(E);
    en = -sum(E.*log(E));
    ens(dim_L+1) = en;
end
%plot entanglement 
figure(3);
n=0:L;
plot(n,ens,'-o');
xlabel('n');
ylabel('entanglement');

