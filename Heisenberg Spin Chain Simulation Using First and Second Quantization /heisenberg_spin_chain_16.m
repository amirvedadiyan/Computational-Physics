clear; 
clc; 
close all; 

%------------initialization------------
L = 14;
J1 = 1; % Nearest neighbor Heisenberg coupling
J2 = 0.3; % Next nearest neighbor Heisenberg coupling
bnd = 1; % boundary condition : 1--> PBC  0 --> OBC  -1 --> APBC
hbar = 1;

%------------pauli matrices------------
sigma_0 = speye(2);
sigma_z = sparse([1,0; 0,-1]);
sigma_x = sparse([0, 1; 1,0]);
sigma_y = sparse([0,-1i; 1i,0]);
sigma_p = (sigma_x + 1i * sigma_y)/2;
sigma_m = (sigma_x - 1i * sigma_y)/2;
n = (sigma_p * sigma_m); % number operator

%------------general problem-> spin operators for n sites------------
sx{1} = hbar/2 * sigma_x;
sy{1} = hbar/2 * sigma_y;
sz{1} = hbar/2 * sigma_z;
ni{1} = hbar/2 * n;
id = sigma_0;
for num_sites = 2:L
    for pos = 1:num_sites-1
        sx{pos} = kron(sx{pos}, sigma_0);
        sy{pos} = kron(sy{pos}, sigma_0);
        sz{pos} = kron(sz{pos}, sigma_0);
        ni{pos} = kron(ni{pos}, sigma_0);
    end
    sx{num_sites} = kron(id, hbar/2*sigma_x);
    sy{num_sites} = kron(id, hbar/2*sigma_y);
    sz{num_sites} = kron(id, hbar/2*sigma_z);
    ni{num_sites} = kron(id,hbar/2*n);
    id = kron(id, sigma_0);
end

%------------Hamiltonian------------
H = sparse(0);
% nearest neighbor interaction
for pos = 1:L-1
    H = H + J1 * (sx{pos}*sx{pos+1} + sy{pos}*sy{pos+1} + sz{pos}*sz{pos+1});
end
if bnd == 1
    H = H + bnd * J1 * (sx{1}*sx{L} + sy{1}*sy{L} + sz{1}*sz{L});
end
% next nearest neighbor interaction
for pos = 1:L-2
    H = H + J2 * (sx{pos}*sx{pos+2} + sy{pos}*sy{pos+2} + sz{pos}*sz{pos+2});
end
if bnd == 1
    H = H + bnd * J2 * (sx{1}*sx{L-1} + sy{1}*sy{L-1} + sz{1}*sz{L-1} + sx{2}*sx{L} + sy{2}*sy{L} + sz{2}*sz{L});
end

%------------U(1) symmtery & total number operator------------
N_tot = sparse(0);
for i =1:L
    N_tot = N_tot + ni{i};
end
N_tot_diag = diag(N_tot);
%States with Sz{tot} = 0 is N_tot = L/2*h/2
F = N_tot_diag==(L*hbar/4);
H_tilde = H(F,F);
[eigenvectors,eigenvalues] = eigs(H_tilde,size(H_tilde,1),'sa');
E = diag(eigenvalues);

%plot lowest 20 energy
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

%------------correlaiton function------------
u = sparse([], [], [], 2^L, 2^L);
u(F,F) = eigenvectors;
[row_indices, ~] = find(F);
% Extract the first row index
first_psi_index = min(row_indices);
psi = u(:,first_psi_index);

SS_corr = zeros(L);
for i1 = 1:L
    for i2 = i1+0:L 
        O = sx{i1}*sx{i2}+sy{i1}*sy{i2}+sz{i1}*sz{i2};
        SS_corr(i1,i2) = psi' * O * psi; 
        SS_corr(i2,i1) = SS_corr(i1,i2);
    end
end

%plot correlation
figure(2);
x = (1:L)-round(L/2);
y = SS_corr(round(L/2),:);
plot(x,y,'-d');
xlabel('distance')
ylabel('Spin-Spin correlation')


%------------entanglement entropy------------
ens = zeros(L,1);
for dim_L = 0:L
    dim_R = L-dim_L;
    psi_tilde = reshape(psi,2^dim_R,2^dim_L);
    rho_L = conj(psi_tilde'*psi_tilde);
    [V,E] = eigs(rho_L,size(rho_L,1),'sa');
    E(E < 1e-13) = 1e-13;
    E = diag(E);
    en = -sum(E.*log(E));
    ens(dim_L) = en;
end
%plot entanglement 
figure(3);
n=0:L;
plot(n,ens,'-o');
xlabel('n');
ylabel('entanglement');
