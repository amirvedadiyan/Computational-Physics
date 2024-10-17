clear; 
clc; 
close all; 

%-------initialization-------
L = 16;
J1 = 1; % Nearest neighbor Heisenberg coupling
bnd = 1; % boundary condition : 1--> PBC  0 --> OBC  -1 --> APBC
hbar = 1;

%-------pauli matrices-------
sigma_0 = speye(2);
sigma_z = sparse([1,0; 0,-1]);
sigma_x = sparse([0, 1; 1,0]);
sigma_y = sparse([0,-1i; 1i,0]);

%-------general problem-> spin operators for n sites-------
sx{1} = hbar/2 * sigma_x;
sy{1} = hbar/2 * sigma_y;
sz{1} = hbar/2 * sigma_z;
id = sigma_0;
for num_sites = 2:L
    for pos = 1:num_sites-1
        sx{pos} = kron(sx{pos}, sigma_0);
        sy{pos} = kron(sy{pos}, sigma_0);
        sz{pos} = kron(sz{pos}, sigma_0);
    end
    sx{num_sites} = kron(id, hbar/2*sigma_x);
    sy{num_sites} = kron(id, hbar/2*sigma_y);
    sz{num_sites} = kron(id, hbar/2*sigma_z);
    id = kron(id, sigma_0);
end
%-------Hamiltonian-------
H = sparse(0);
% nearest neighbor interaction only between (1,2),(3,4),(5,6)...
for pos = 1:2:L
    H = H + J1 * (sx{pos}*sx{pos+1}+sy{pos}*sy{pos+1}+sz{pos}*sz{pos+1});
end
%-------diagonalization-------
n_eigs = min(10, size(H,1));
[U, D] = eigs(H, n_eigs, 'sa');
E = diag(D);
psi_0 = U(:,1); % |psi> ground-state wavefunction

%-------time evolution-------
start_t = 0;
end_t = 10;
step = 0.05;
time = start_t:step:end_t;
%--------New Hamiltonian-----
H_t = sparse(0);
for pos = 1:L-1
    H_t = H_t + J1 * (sx{pos}*sx{pos+1} + sy{pos}*sy{pos+1} + sz{pos}*sz{pos+1});
end
if bnd == 1
    H_t = H_t + bnd * J1 * (sx{1}*sx{L} + sy{1}*sy{L} + sz{1}*sz{L});
end
psi_t = psi_0;
ee = zeros(length(time),1);
SS_corr4 = {};
SS_corr4p = {};
mi=1;
si=1;
for t = 1:length(time)
    psi_t = (speye(2^L) - 1i * hbar * step * H_t) * psi_t;
    %normalization
    norm_psi = norm(psi_t);
    % Normalize the sparse array
    psi_t = psi_t / norm_psi;
    % norm_psi_t = norm(psi_t); it will be 1
    % calculate correlation function <Si(t).Sj(t)> & <Si(t).Sj(0)> 
    if t == 1 || t == 2/step || t == 4/step || t == 6/step
        SS_corr1 = zeros(L);
        SS_corr2 = zeros(L);
        for i1 = 1:L
            for i2 = i1+0:L 
                O = sx{i1}*sx{i2}+sy{i1}*sy{i2}+sz{i1}*sz{i2};
                SS_corr1(i1,i2) = psi_t' * O * psi_t; 
                SS_corr1(i2,i1) = SS_corr1(i1,i2);
            end
        end
        SS_corr4{si} = SS_corr1;
        for i1 = 1:L
            for i2 = i1+0:L 
                O = sx{i1}*sx{i2}+sy{i1}*sy{i2}+sz{i1}*sz{i2};
                SS_corr2(i1,i2) = psi_t' * O * psi_0 ; 
                SS_corr2(i2,i1) = SS_corr2(i1,i2);
            end
        end
        SS_corr4p{si} = SS_corr2;
        si = si + 1;
    end
    %-------calculate entanglement-------
    dim_R = L/2;
    dim_L = L/2;
    psi_tilde = reshape(psi_t,2^dim_R,2^dim_L);
    rho_L = conj(psi_tilde'*psi_tilde);
    [V,E] = eigs(rho_L,size(rho_L,1),'sa');
    E(E < 1e-13) = 1e-13;
    E = diag(E);
    ee(mi) = -sum(E.*log(E));
    mi = mi +1;
    
    
end

% plot correlation functions <Si(t).Sj(t)> in times: 2s,5s,8s
figure(1);
x = (1:L)-round(L/2);
y1 = SS_corr4{1};
y2 = SS_corr4{2};
y3 = SS_corr4{3};
y4 = SS_corr4{4};
y1 = y1(round(L/2),:);
y2 = y2(round(L/2),:);
y3 = y3(round(L/2),:);
y4 = y4(round(L/2),:);
plot(x,y1,'-d');
hold on;
plot(x,y2,'-d');
plot(x,y3,'-d');
plot(x,y4,'-d');
xlabel('distance')
ylabel('Spin-Spin correlation <Si(t).Sj(t)>')
legend('t=0','t=2' ,'t=4','t=6')
% plot correlation <Si(t).Sj(0)> in times: 2s,5s,8s
figure(2);
y1p = SS_corr4p{1};
y2p = SS_corr4p{2};
y3p = SS_corr4p{3};
y4p = SS_corr4p{4};
y1p = y1p(round(L/2),:);
y2p = y2p(round(L/2),:);
y3p = y3p(round(L/2),:);
y4p = y4p(round(L/2),:);
plot(x,y1,'-d');
hold on;
plot(x,y2p,'-d');
plot(x,y3p,'-d');
plot(x,y4p,'-d');
xlabel('distance')
ylabel('Spin-Spin correlation <Si(t).Sj(0)>')
legend('t=0','t=2' ,'t=4','t=6')


% plot entanglement with times
figure(3)
plot(time,ee,'-o');
xlabel('time')
ylabel('Entanglement')

