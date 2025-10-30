%{
Copyright (c) 2025 Qiang Gao
Author: Qiang Gao <gq201277@gmail.com>
Description: Bootstrap the stiffness of a flatband superconductor
%}

function [Stiffness_boot,runtime] = bootstrapSC(Nx, Ny, N, A)
checker(Nx,Ny,N,A)
format long
tic
kx0 = -1/2+1/Nx/2:1/Nx:1/2;
ky0 = -1/2+1/Ny/2:1/Ny:1/2;

a1 = [1, 0];
a2 = [0, 1];

A_mat = [a1; a2];
B_mat = 2*pi*inv(A_mat)';

[KX, KY] = meshgrid(kx0, ky0);
k0_ind   = [KX(:), KY(:)];
k0 = k0_ind*B_mat;

Nk = length(k0_ind);

[M_G2D,M_rho,DkkBZ,D2_mask_total,D2_mask,M_k_Qmk,M_T2_to_D_selected,M_T2_to_rho_selected,Extra_part_T2_selected]=build_cons(k0_ind,N);

[F_Q_kkp_ssp_A,F_k_sz_A] = H_TFB(Nk,A,k0,DkkBZ);

num_Sz_T2 = 10;
num_Sz_D = 3;
num_Sz_Q = num_Sz_D;
num_Sz_G = 6;
num_Sz_dc = 2;

ops = sdpsettings('solver','mosek','verbose',0);
M_D_s = sdpvar(Nk,Nk,Nk,2,'hermitian','complex');
M_D_s = D2_mask.*M_D_s;

M_D_s_conj = M_D_s(:,:,:,1);

for Qi = 1:Nk
    M_D_s_conj(:,:,Qi) = conj(M_D_s(:,:,Qi,1));
end

M_D = cat(4,M_D_s(:,:,:,1),M_D_s_conj,M_D_s(:,:,:,2));
rho_k_sz = sdpvar(Nk,num_Sz_dc,'full','real');

Constraints = [rho_k_sz == reshape(M_rho*reshape(M_D,[],1),[],2),sum(sum(rho_k_sz)) == N];

for sz = 1:num_Sz_dc
    for i = 1:Nk
        Constraints = [Constraints, rho_k_sz(i,sz) >= 0];
    end
end

M_G = reshape(M_G2D*reshape(M_D,[],1),Nk,Nk,Nk,num_Sz_G);

M_Q = M_D;
Sz_D_map = [1,1;2,2;1,2];
for Sz = 1:num_Sz_Q
    for Qi = 1:Nk
        M_Q(:,:,Qi,Sz) = M_D(:,:,Qi,Sz) ...
            +D2_mask_total(:,:,Qi,Sz).*diag(1-rho_k_sz(:,Sz_D_map(Sz,1))-M_k_Qmk(:,:,Qi)*rho_k_sz(:,Sz_D_map(Sz,2)));
    end
end

M_T2 = reshape(M_T2_to_D_selected*reshape(M_D,[],1)+M_T2_to_rho_selected*reshape(rho_k_sz,[],1)...
    +Extra_part_T2_selected,[Nk^2,Nk^2,1,num_Sz_T2]);



% D condition
for Sz = 1:num_Sz_D
    for Qi = 1:Nk
        Constraints = [Constraints, M_D(:,:,Qi,Sz)>=0];
    end
end

% G condition
for Sz = 1:2
    for Qi = 1:Nk
        Constraints = [Constraints, (M_G(:,:,Qi,Sz)+M_G(:,:,Qi,Sz)')/2 + diag(rho_k_sz(:,mod(Sz-1,2)+1))>=0];
    end
end

for Qi = 1:Nk
    Constraints = [Constraints, ...
        [(M_G(:,:,Qi,3)+M_G(:,:,Qi,3)')/2 + diag(rho_k_sz(:,1)),M_G(:,:,Qi,5);...
        M_G(:,:,Qi,6),(M_G(:,:,Qi,4)+M_G(:,:,Qi,4)')/2 + diag(rho_k_sz(:,2))] >=0];
end

% Q condition
M_Q = M_D;
for Sz = 1:num_Sz_Q
    for Qi = 1:Nk
        Constraints = [Constraints, (M_Q(:,:,Qi,Sz)+M_Q(:,:,Qi,Sz)')/2 >=0];
    end
end

% T2 condition
Constraints = [Constraints,M_T2(:,:,1,1)+M_T2(:,:,1,1)'>=0,...
    [M_T2(:,:,1,2),M_T2(:,:,1,3);M_T2(:,:,1,4),M_T2(:,:,1,5)]>=0];


Objective = real(transpose(F_Q_kkp_ssp_A)*reshape(M_G(:,:,:,3:6),[],1)+ reshape(F_k_sz_A,1,[])*reshape(rho_k_sz,[],1));

diagnostic = optimize(Constraints, Objective, ops);
if diagnostic.problem ~= 0
    fprintf('The solver did not give a converged result on this data point: Nx=%d, Ny=%d, N=%d',Nx,Ny,N)
end
M_G_value_T2 = value(M_G);
rho_k_sz_value_T2 = value(rho_k_sz);
Energy = transpose(F_Q_kkp_ssp_A)*reshape(M_G_value_T2(:,:,:,3:6),[],1)+ reshape(F_k_sz_A,1,[])*reshape(rho_k_sz_value_T2,[],1);
yalmip('clear')
runtime = toc;
Stiffness_boot = Energy/(A(1)^2+A(2)^2)/(Nx*Ny)^2/4/4;
assert(abs(imag(Stiffness_boot))<1e-10)
Stiffness_boot = real(Stiffness_boot);

end

function checker(Nx,Ny,N,A)
if ~(isfinite(Nx) && (Nx == fix(Nx))) || ~(isfinite(Ny) && (Ny == fix(Ny)))  || Nx<=0 || Ny <=0
    error('Invalid inputs: Nx and Ny must be positive integers')
elseif mod(Nx,2)~=1 || mod(Ny,2)~=1 
    error('This algorithm only works for odd by odd lattices')
end
if ~(isfinite(N) && (N == fix(N))) || N<=0
    error('Invalid inputs: the particle number N must be a positive integer')
elseif mod(N,2)~=0
    error('Please enter an even number for N')
end
if norm(A)>0.1
    warning('Please choose a smaller gauge insertion for this algorithm to work better.')
end
end

function [F_Q_kkp_ssp_A,F_k_sz_A] = H_TFB(Nk,A,k0,DkkBZ)

D_AB_k = zeros(Nk,4,2,2);

D_AB_k(:,:,:,1) = ones(Nk,4,2);

for q = 1:Nk
    for pm = 1:2
        qx = k0(q,1);
        qy = k0(q,2);
        D_AB_k(q,1,pm,2) = (-1)^(pm-1);
        D_AB_k(q,2,pm,2) = (-1)^(pm-1)*exp(-1i*qx);
        D_AB_k(q,3,pm,2) = (-1)^(pm-1)*exp(-1i*qy);
        D_AB_k(q,4,pm,2) = (-1)^(pm-1)*exp(-1i*qx-1i*qy);
    end
end

H_k_single_pA = zeros(2,2,Nk);
H_k_single_mA = zeros(2,2,Nk);
t1 = 1;
t2 = 1/sqrt(2);
t5 = (1-sqrt(2))/4;

for i = 1:Nk
    kx = k0(i,1)+A(1);
    ky = k0(i,2)+A(2);
    B_k_0 = -2*t5*(cos(2*kx)+cos(2*ky));
    B_k_x_iy = -t1*(exp(1i*(kx+ky+pi/4))+exp(1i*pi/4)+exp(1i*(kx-pi/4))+exp(1i*(ky-pi/4)));
    B_k_z = -2*t2*(cos(ky)-cos(kx));

    H_k_single_pA(:,:,i) = [B_k_0+B_k_z, B_k_x_iy; conj(B_k_x_iy), B_k_0-B_k_z];
end

for i = 1:Nk
    kx = k0(i,1)-A(1);
    ky = k0(i,2)-A(2);

    B_k_0 = -2*t5*(cos(2*kx)+cos(2*ky));
    B_k_x_iy = -t1*(exp(1i*(kx+ky+pi/4))+exp(1i*pi/4)+exp(1i*(kx-pi/4))+exp(1i*(ky-pi/4)));
    B_k_z = -2*t2*(cos(ky)-cos(kx));

    H_k_single_mA(:,:,i) = [B_k_0+B_k_z, B_k_x_iy; conj(B_k_x_iy), B_k_0-B_k_z];
end

energy_k_single_pA = zeros(2,Nk);
energy_k_single_mA = zeros(2,Nk);

U_k_alpha_pA = zeros(2,2,Nk);
U_k_alpha_mA = zeros(2,2,Nk);

for i = 1:Nk
    [V,D] = eig(H_k_single_pA(:,:,i));
    [energy_sorted,index_sort_energy] = sort(real(diag(D)));
    energy_k_single_pA(:,i) = energy_sorted;
    U_k_alpha_pA(:,:,i) = V(:,index_sort_energy);
end

for i = 1:Nk
    [V,D] = eig(H_k_single_mA(:,:,i));
    [energy_sorted,index_sort_energy] = sort(real(diag(D)));
    energy_k_single_mA(:,i) = energy_sorted;
    U_k_alpha_mA(:,:,i) = V(:,index_sort_energy);
end


U_k_alpha_flat_pA = reshape(U_k_alpha_pA(:,1,:),[],Nk);
U_k_alpha_flat_mA = reshape(U_k_alpha_mA(:,1,:),[],Nk);

for i = 1:Nk
    U_k_alpha_flat_pA(:,i) = U_k_alpha_flat_pA(:,i)*conj(U_k_alpha_flat_pA(2,i))/abs(conj(U_k_alpha_flat_pA(2,i)));
    U_k_alpha_flat_mA(:,i) = U_k_alpha_flat_mA(:,i)*conj(U_k_alpha_flat_mA(2,i))/abs(conj(U_k_alpha_flat_mA(2,i)));
end

U_k_alpha_flat_spin_up = U_k_alpha_flat_pA;
U_k_alpha_flat_spin_down = conj(flip(U_k_alpha_flat_mA,2));

U_k_alpha_flat_spin = cat(3,U_k_alpha_flat_spin_up,U_k_alpha_flat_spin_down);

num_Sz_F = 4;
Sz_F_map = [1,1;2,2;1,2;2,1]; 
k_index_ordered_F = zeros(Nk^3*num_Sz_F,4);

count = 0;
for Sz = 1:num_Sz_F
    for Qi = 1:Nk
        for kp = 1:Nk
            for k = 1:Nk
                count = count + 1;
                k_index_ordered_F(count,:) = [k,kp,Qi,Sz+2];
            end
        end
    end
end


F_Q_kkp_ssp_A = zeros(length(k_index_ordered_F),1);
F_k_sz_A = zeros(Nk,2);

for i = 1:length(F_Q_kkp_ssp_A)
    k1 = k_index_ordered_F(i,1);
    k1p = k_index_ordered_F(i,2);
    Qi = k_index_ordered_F(i,3);
    Sz = k_index_ordered_F(i,4);
    sz = Sz_F_map(Sz-2,1);
    szp = Sz_F_map(Sz-2,2);

    xi_sz = (-1)^(sz-1);
    xi_szp = (-1)^(szp-1);

    k2 = DkkBZ(k1,Qi);
    k2p = DkkBZ(k1p,Qi);

    for alpha = 1:2
        for beta = 1:2
            for pm = 1:2
                for I = 1:4
                    F_Q_kkp_ssp_A(i) = F_Q_kkp_ssp_A(i) + 1/2*xi_sz*xi_szp...
                        *D_AB_k(DkkBZ(k1,k2),I,pm,alpha)*D_AB_k(DkkBZ(k2p,k1p),I,pm,beta)...
                        *(conj(U_k_alpha_flat_spin(alpha,k1,sz))*U_k_alpha_flat_spin(alpha,k2,sz)...
                        *conj(U_k_alpha_flat_spin(beta,k2p,szp))*U_k_alpha_flat_spin(beta,k1p,szp));
                end
            end
        end
    end

end

for sz = 1:2
    for i = 1:Nk
        k1 = i;
        k1p = i;
        szp = sz;

        temp = 0;
        for Qi = 1:Nk
            k2 = DkkBZ(k1,Qi);
            k2p = k2;
            for alpha = 1:2
                for beta = 1:2
                    for pm = 1:2
                        for I = 1:4
                            temp = temp + 1/2 ...
                                *D_AB_k(DkkBZ(k1,k2),I,pm,alpha)*D_AB_k(DkkBZ(k2p,k1p),I,pm,beta)...
                                *(conj(U_k_alpha_flat_spin(alpha,k1,sz))*U_k_alpha_flat_spin(alpha,k2,sz)...
                                *conj(U_k_alpha_flat_spin(beta,k2p,szp))*U_k_alpha_flat_spin(beta,k1p,szp));
                        end
                    end
                end
            end
        end
        F_k_sz_A(i,sz) = temp;
    end
end

end

function [M_G2D,M_d_dc_c,DkkBZ,D2_mask_total,D2_mask_SU2_total,M_k_Qmk,M_T2_to_D_selected,M_T2_to_rho_selected,Extra_part_T2_selected]=build_cons(k0_ind,N)
Nk = length(k0_ind);
[M_G2D,M_d_dc_c,SkkBZ,DkkBZ,D2_mask_total,D2_mask_SU2_total,M_k_Qmk] = construct_DQG(k0_ind,N);
[M_T2_to_D_selected,M_T2_to_rho_selected,Extra_part_T2_selected] = construct_T2(Nk,SkkBZ,DkkBZ);
end

function [M_g2d,M_d2rho,SkkBZ,DkkBZ,D2_mask_total,D2_mask_SU2_total,M_k_Qmk] = construct_DQG(k0_ind,N)

[SkkBZ,DkkBZ,D2_mask_total,D2_mask_SU2_total,M_k_Qmk]=lattice_mat(k0_ind);
Nk = length(k0_ind);
M_g2d = G2D(Nk,SkkBZ,DkkBZ);
M_d2rho = D2rho(Nk,N,DkkBZ);
end

function [index,index_finder] = D_indexing(Nk)
num_Sz = 3;
k_index_ordered = zeros(Nk^3*num_Sz,4);
index_finder = zeros(Nk,Nk,Nk,num_Sz);

count = 0;
for Sz = 1:num_Sz
    for Qi = 1:Nk
        for kp = 1:Nk
            for k = 1:Nk
                count = count + 1;
                k_index_ordered(count,:) = [k,kp,Qi,Sz];
                index_finder(k,kp,Qi,Sz) = count;
            end
        end
    end
end
index = k_index_ordered;
end

function [index_g,index_finder_g] = G_indexing(Nk)
num_Sz = 6;
k_index_ordered = zeros(Nk^3*num_Sz,4);
index_finder = zeros(Nk,Nk,Nk,num_Sz);

count = 0;
for Sz = 1:num_Sz
    for Qi = 1:Nk
        for kp = 1:Nk
            for k = 1:Nk
                count = count + 1;
                k_index_ordered(count,:) = [k,kp,Qi,Sz];
                index_finder(k,kp,Qi,Sz) = count;
            end
        end
    end
end

index_g = k_index_ordered;

if nargout>=2
    index_finder_g = index_finder;
end
end

function [Skk,Dkk,D2_mask_total,D2_mask_SU2_total,M_k_Qmk]=lattice_mat(k0_ind)
Nk = length(k0_ind);
Skk = zeros(Nk,Nk);
for i = 1:Nk
    for j = 1:Nk
        temp1 = mod(k0_ind(i,1)+k0_ind(j,1),1);
        if temp1 > 1/2
            temp1 = temp1 - 1;
        end
        temp2 = mod(k0_ind(i,2)+k0_ind(j,2),1);
        if temp2 > 1/2
            temp2 = temp2 - 1;
        end
        for m = 1:Nk
            if norm([temp1,temp2]-k0_ind(m,:))<0.0001
                Skk(i,j) = m;
            end
        end
    end
end
Skk = round(Skk);

Dkk = zeros(Nk,Nk);
for i = 1:Nk
    for j = 1:Nk
        temp1 = mod(k0_ind(i,1)-k0_ind(j,1),1);
        if temp1 > 1/2
            temp1 = temp1 - 1;
        end
        temp2 = mod(k0_ind(i,2)-k0_ind(j,2),1);
        if temp2 > 1/2
            temp2 = temp2 - 1;
        end
        for m = 1:Nk
            if norm([temp1,temp2]-k0_ind(m,:))<0.0001
                Dkk(i,j) = m;
            end
        end
    end
end
D2_uuuu_mask = ones(Nk,Nk,Nk);

for Qi = 1:Nk
    for k = 1:Nk
        for kp = 1:Nk
            if Dkk(Qi,k) <= k || Dkk(Qi,kp) <= kp
                D2_uuuu_mask(k,kp,Qi) = 0;
            end
        end
    end
end

D2_dddd_mask = D2_uuuu_mask;
D2_uddu_mask = ones(Nk,Nk,Nk);
D2_mask_total = cat(4,D2_uuuu_mask,D2_dddd_mask,D2_uddu_mask);
D2_mask_SU2_total = cat(4,D2_uuuu_mask,D2_uddu_mask);

M_k_Qmk = zeros(Nk,Nk,Nk);
for Qi = 1:Nk
    for k = 1:Nk
        M_k_Qmk(k,Dkk(Qi,k),Qi) = 1;
    end
end
end

function rho_k_sz = D2rho(Nk,N,DkkBZ)
    [index_ddcc,index_finder_ddcc] = D_indexing(Nk);
    dim_ddcc = length(index_ddcc);

    num_sz = 2;
    ksz_index = zeros(num_sz*Nk,2); 
    count = 0;
    for sz = 0:num_sz-1
        for i = 1:Nk
            count = count + 1;
            ksz_index(count,:) = [i,sz];
        end
    end

    dim_rho = length(ksz_index);

    coords = zeros(6*Nk*dim_rho,2);
    values = zeros(6*Nk*dim_rho,1);

    count = 0;

    for i = 1:dim_rho
        k = ksz_index(i,1);
        sz = ksz_index(i,2);

        for szt = 0:num_sz-1
            for Qt = 1:Nk
                k1t = k;
                k2t = DkkBZ(Qt,k);
                k1tp = k;
                k2tp = DkkBZ(Qt,k);

                    sz1t = sz;
                    sz2t = szt;
                    sz1tp = sz;
                    sz2tp = szt;

                    if sz==szt

                        if k1t ~= k2t && k1tp ~= k2tp
                        Sz_ddcc = sz+1;

                        parity = 1;
                        parity_p = 1;

                        if k1t > k2t
                            k1t = k2t;
                            parity = -1;
                        end

                        if k1tp > k2tp
                            k1tp = k2tp;
                            parity_p = -1;
                        end

                        count = count + 1;
                        coords(count,:) = [i,index_finder_ddcc(k1t,k1tp,Qt,Sz_ddcc)];
                        values(count) = parity*parity_p/(N-1);
                        end
                    else
                        Sz_ddcc = 3;
                        parity = 1;
                        parity_p = 1;
                        if sz1t == 1
                            assert(sz2t == 0)
                            k1t = k2t;
                            parity = -1;
                        end

                        if sz1tp == 1
                            assert(sz2tp == 0)
                            k1tp = k2tp;
                            parity_p = -1;
                        end

                        count = count + 1;
                        coords(count,:) = [i,index_finder_ddcc(k1t,k1tp,Qt,Sz_ddcc)];
                        values(count) = parity*parity_p/(N-1);
                    end
            end
        end

    end

    coords = coords(abs(values)>0,:);
    values = values(abs(values)>0);
    rho_k_sz = sparse(coords(:,1),coords(:,2),values,dim_rho,dim_ddcc);


end

function M_g2d = G2D(Nk,SkkBZ,DkkBZ)
[index_d,index_finder_d] = D_indexing(Nk);
index_g = G_indexing(Nk);

dim_d = length(index_d);
dim_g = length(index_g);

Sz_G = [0,1,1,0; 1,0,0,1; 0,0,0,0; 1,1,1,1; 0,0,1,1; 1,1,0,0];

coords = zeros(dim_g,2);
values = zeros(dim_g,1);

count = 0;

for i = 1:dim_g
    k1 = index_g(i,1);
    k1p = index_g(i,2);
    Qi = index_g(i,3);

    k2 = DkkBZ(k1,Qi);
    k2p = DkkBZ(k1p,Qi);

    Sz_g = index_g(i,4);

    sz1 = Sz_G(Sz_g,1);
    sz1p = Sz_G(Sz_g,4);
    sz2 = Sz_G(Sz_g,2);
    sz2p = Sz_G(Sz_g,3);

    k1t = k1;
    k2t = k2p;
    k1tp = k1p;
    k2tp = k2;

    sz1t = sz1;
    sz2t = sz2p;
    sz1tp = sz1p;
    sz2tp = sz2;

    assert(sz1t+sz2t == sz1tp + sz2tp)
    assert(SkkBZ(k2t,k1t) == SkkBZ(k2tp,k1tp))
    Qi_d = SkkBZ(k2t,k1t);
    if Sz_g == 3 || Sz_g == 4
        if k2t ~= k1t && k2tp ~= k1tp
            Sz_d = Sz_g-2;

            parity = 1;
            parity_p = 1;

            if k1t > k2t
                k1t = k2t;
                parity = -1;
            end

            if k1tp > k2tp
                k1tp = k2tp;
                parity_p = -1;
            end

            count = count + 1;
            coords(count,:) = [i,index_finder_d(k1t,k1tp,Qi_d,Sz_d)];
            values(count) = -parity*parity_p;
        end

    else
        Sz_d = 3;

        parity = 1;
        parity_p = 1;
        if sz1t == 1
            assert(sz2t == 0)
            k1t = k2t;
            parity = -1;
        end

        if sz1tp == 1
            assert(sz2tp == 0)
            k1tp = k2tp;
            parity_p = -1;
        end

        count = count + 1;
        coords(count,:) = [i,index_finder_d(k1t,k1tp,Qi_d,Sz_d)];
        values(count) = -parity*parity_p;
    end

end

coords = coords(abs(values)>0,:);
values = values(abs(values)>0);

M_g2d = sparse(coords(:,1),coords(:,2),values,dim_g,dim_d);
end

function [M_T2_to_D,M_T2_to_rho,Extra_part_T2] = construct_T2(Nk,SkkBZ,DkkBZ)
[index_T1,~,index_finder_T1] = T1_index(Nk);
[index_d,index_finder_d] = D_indexing(Nk);
[index_rho,index_finder_rho] = rho_index(Nk);
num_Sz_T2 = 10;
Q_T2 = (Nk+1)/2;
k1_k2_ordered_finder = zeros(Nk,Nk);

count = 0;
for k1 = 1:Nk-1
    for k2 = k1+1:Nk
        count = count + 1;
        k1_k2_ordered_finder(k1,k2) = count;
    end
end
assert(count==Nk*(Nk-1)/2)

T2_k1_k2 = zeros(Nk^2,2);

count = 0;
for k1 = 1:Nk
    for k2 = 1:Nk
        count = count + 1;
        T2_k1_k2(count,:) = [k1,k2];
    end
end

dim_T2 = Nk^2*Nk^2*length(Q_T2)*num_Sz_T2;
dim_T1 = length(index_T1);
dim_d = length(index_d);
dim_rho = length(index_rho);

coords_to_T1 = zeros(dim_T2,2);
values_to_T1 = zeros(dim_T2,1);
coords_to_D = zeros(dim_T2,2);
values_to_D = zeros(dim_T2,1);
coords_to_rho = zeros(dim_T2,2);
values_to_rho = zeros(dim_T2,1);

count_T1 = 0;
count_D = 0;
count_rho = 0;


Extra_part_T2 = zeros(dim_T2,1);



i = 0;
for Sz_T2 = 1:num_Sz_T2
    for Qi = 1:length(Q_T2)
        for kp = 1:Nk^2
            for k = 1:Nk^2
                i = i+1;
                Qi_T2 = Q_T2(Qi);
                k1_ind = T2_k1_k2(k,1);
                k2_ind = T2_k1_k2(k,2);
                k3_ind = DkkBZ(SkkBZ(k1_ind,k2_ind),Qi_T2);
                k1p_ind = T2_k1_k2(kp,1);
                k2p_ind = T2_k1_k2(kp,2);
                k3p_ind = DkkBZ(SkkBZ(k1p_ind,k2p_ind),Qi_T2);

                if Sz_T2 == 1
                    Sz_T1 = 3;
                    if k2_ind > k1_ind && k2p_ind > k1p_ind
                        Qi_T1 = SkkBZ(k1_ind,SkkBZ(k2_ind,k3p_ind));
                        count_T1 = count_T1 + 1;
                        coords_to_T1(count_T1,:) = [i,index_finder_T1(k1_k2_ordered_finder(k1_ind,k2_ind),k1_k2_ordered_finder(k1p_ind,k2p_ind),Qi_T1,Sz_T1)];
                        values_to_T1(count_T1) = -1;

                        if k3_ind == k3p_ind
                            Qi_d = SkkBZ(k1_ind,k2_ind);
                            Sz_d = 1;
                            sz_rho = 0;

                            count_D = count_D + 1;
                            coords_to_D(count_D,:) = [i,index_finder_d(k1_ind,k1p_ind,Qi_d,Sz_d)];
                            values_to_D(count_D) = 2;

                            if k1_ind == k1p_ind && k2_ind == k2p_ind
                                Extra_part_T2(i) = Extra_part_T2(i) + 1;
                                count_rho = count_rho+1;
                                coords_to_rho(count_rho,:) = [i,index_finder_rho(k1_ind,sz_rho+1)];
                                values_to_rho(count_rho) = -1;

                                count_rho = count_rho+1;
                                coords_to_rho(count_rho,:) = [i,index_finder_rho(k2_ind,sz_rho+1)];
                                values_to_rho(count_rho) = -1;
                            end

                        end
                    end

                end

                if Sz_T2 == 10
                    Sz_T1 = 4;
                    if k2_ind > k1_ind && k2p_ind > k1p_ind
                        Qi_T1 = SkkBZ(k1_ind,SkkBZ(k2_ind,k3p_ind));
                        count_T1 = count_T1 + 1;
                        coords_to_T1(count_T1,:) = [i,index_finder_T1(k1_k2_ordered_finder(k1_ind,k2_ind),k1_k2_ordered_finder(k1p_ind,k2p_ind),Qi_T1,Sz_T1)];
                        values_to_T1(count_T1) = -1;

                        if k3_ind == k3p_ind
                            Qi_d = SkkBZ(k1_ind,k2_ind);
                            Sz_d = 2;
                            sz_rho = 1;

                            count_D = count_D + 1;
                            coords_to_D(count_D,:) = [i,index_finder_d(k1_ind,k1p_ind,Qi_d,Sz_d)];
                            values_to_D(count_D) = 2;

                            if k1_ind == k1p_ind && k2_ind == k2p_ind
                                Extra_part_T2(i) = Extra_part_T2(i) + 1;
                                count_rho = count_rho+1;
                                coords_to_rho(count_rho,:) = [i,index_finder_rho(k1_ind,sz_rho+1)];
                                values_to_rho(count_rho) = -1;

                                count_rho = count_rho+1;
                                coords_to_rho(count_rho,:) = [i,index_finder_rho(k2_ind,sz_rho+1)];
                                values_to_rho(count_rho) = -1;
                            end

                        end
                    end

                end

                if Sz_T2 == 2
                    Sz_T1 = 1;
                    if k2_ind > k1_ind && k2p_ind > k1p_ind
                        if k3_ind == k3p_ind
                            Qi_d = SkkBZ(k1_ind,k2_ind);
                            Sz_d = 1;
                            sz_rho = 0;

                            count_D = count_D + 1;
                            coords_to_D(count_D,:) = [i,index_finder_d(k1_ind,k1p_ind,Qi_d,Sz_d)];
                            values_to_D(count_D) = 2;

                            if k1_ind == k1p_ind && k2_ind == k2p_ind
                                Extra_part_T2(i) = Extra_part_T2(i) + 1;
                                count_rho = count_rho+1;
                                coords_to_rho(count_rho,:) = [i,index_finder_rho(k1_ind,sz_rho+1)];
                                values_to_rho(count_rho) = -1;

                                count_rho = count_rho+1;
                                coords_to_rho(count_rho,:) = [i,index_finder_rho(k2_ind,sz_rho+1)];
                                values_to_rho(count_rho) = -1;
                            end

                        end

                        if k3p_ind ~= k1_ind && k3p_ind ~= k2_ind && k3_ind ~= k1p_ind && k3_ind ~= k2p_ind
                            Qi_T1 = SkkBZ(k1_ind,SkkBZ(k2_ind,k3p_ind));
                            parity = 1;
                            parity_p = 1;

                            if k3p_ind < k1_ind
                                k2_ind = k1_ind;
                                k1_ind = k3p_ind;
                                parity = 1;
                            elseif k3p_ind < k2_ind
                                k2_ind = k3p_ind;
                                parity = -1;
                            end

                            if k3_ind < k1p_ind
                                k2p_ind = k1p_ind;
                                k1p_ind = k3_ind;
                                parity_p = 1;
                            elseif k3_ind < k2p_ind
                                k2p_ind = k3_ind;
                                parity_p = -1;
                            end

                            count_T1 = count_T1 + 1;
                            coords_to_T1(count_T1,:) = [i,index_finder_T1(k1_k2_ordered_finder(k1_ind,k2_ind),k1_k2_ordered_finder(k1p_ind,k2p_ind),Qi_T1,Sz_T1)];
                            values_to_T1(count_T1) = -parity*parity_p;
                        end
                    end

                end

                if Sz_T2 == 6
                    Sz_T1 = 2;
                    if k2_ind > k1_ind && k2p_ind > k1p_ind
                        if k3_ind == k3p_ind
                            Qi_d = SkkBZ(k1_ind,k2_ind);
                            Sz_d = 2;
                            sz_rho = 1;

                            count_D = count_D + 1;
                            coords_to_D(count_D,:) = [i,index_finder_d(k1_ind,k1p_ind,Qi_d,Sz_d)];
                            values_to_D(count_D) = 2;

                            if k1_ind == k1p_ind && k2_ind == k2p_ind
                                Extra_part_T2(i) = Extra_part_T2(i) + 1;
                                count_rho = count_rho+1;
                                coords_to_rho(count_rho,:) = [i,index_finder_rho(k1_ind,sz_rho+1)];
                                values_to_rho(count_rho) = -1;

                                count_rho = count_rho+1;
                                coords_to_rho(count_rho,:) = [i,index_finder_rho(k2_ind,sz_rho+1)];
                                values_to_rho(count_rho) = -1;
                            end

                        end

                        if k3p_ind ~= k1_ind && k3p_ind ~= k2_ind && k3_ind ~= k1p_ind && k3_ind ~= k2p_ind
                            Qi_T1 = SkkBZ(k1_ind,SkkBZ(k2_ind,k3p_ind));
                            parity = 1;
                            parity_p = 1;

                            if k3p_ind < k1_ind
                                k2_ind = k1_ind;
                                k1_ind = k3p_ind;
                                parity = 1;
                            elseif k3p_ind < k2_ind
                                k2_ind = k3p_ind;
                                parity = -1;
                            end

                            if k3_ind < k1p_ind
                                k2p_ind = k1p_ind;
                                k1p_ind = k3_ind;
                                parity_p = 1;
                            elseif k3_ind < k2p_ind
                                k2p_ind = k3_ind;
                                parity_p = -1;
                            end

                            count_T1 = count_T1 + 1;
                            coords_to_T1(count_T1,:) = [i,index_finder_T1(k1_k2_ordered_finder(k1_ind,k2_ind),k1_k2_ordered_finder(k1p_ind,k2p_ind),Qi_T1,Sz_T1)];
                            values_to_T1(count_T1) = -parity*parity_p;
                        end
                    end

                end

                if Sz_T2 == 5 
                    Sz_T1 = 4; 
                    if k3_ind == k3p_ind
                        Qi_d = SkkBZ(k1_ind,k2_ind);
                        Sz_d = 3;
                        sz_rho_1 = 0;
                        sz_rho_2 = 1;

                        count_D = count_D + 1;
                        coords_to_D(count_D,:) = [i,index_finder_d(k1_ind,k1p_ind,Qi_d,Sz_d)];
                        values_to_D(count_D) = 2;

                        if k1_ind == k1p_ind && k2_ind == k2p_ind
                            Extra_part_T2(i) = Extra_part_T2(i) + 1;
                            count_rho = count_rho+1;
                            coords_to_rho(count_rho,:) = [i,index_finder_rho(k1_ind,sz_rho_1+1)];
                            values_to_rho(count_rho) = -1;

                            count_rho = count_rho+1;
                            coords_to_rho(count_rho,:) = [i,index_finder_rho(k2_ind,sz_rho_2+1)];
                            values_to_rho(count_rho) = -1;
                        end
                    end

                    if k2_ind ~= k3p_ind && k2p_ind ~= k3_ind
                        Qi_T1 = SkkBZ(k1_ind,SkkBZ(k2_ind,k3p_ind));
                        parity = 1;
                        parity_p = 1;

                        if k2_ind > k3p_ind
                            temp = k2_ind;
                            k2_ind = k3p_ind;
                            k3p_ind = temp;
                            parity = -1;
                        end

                        if k2p_ind > k3_ind
                            temp = k2p_ind;
                            k2p_ind = k3_ind;
                            k3_ind = temp;
                            parity_p = -1;
                        end

                        count_T1 = count_T1 + 1;
                        coords_to_T1(count_T1,:) = [i,index_finder_T1(k1_k2_ordered_finder(k2_ind,k3p_ind),k1_k2_ordered_finder(k2p_ind,k3_ind),Qi_T1,Sz_T1)];
                        values_to_T1(count_T1) = -parity*parity_p;

                    end

                end

                if Sz_T2 == 9
                    Sz_T1 = 3; 
                    if k3_ind == k3p_ind
                        Qi_d = SkkBZ(k1_ind,k2_ind);
                        Sz_d = 3;
                        sz_rho_2 = 0; 
                        sz_rho_1 = 1;

                        count_D = count_D + 1;
                        coords_to_D(count_D,:) = [i,index_finder_d(k2_ind,k2p_ind,Qi_d,Sz_d)];
                        values_to_D(count_D) = 2;

                        if k1_ind == k1p_ind && k2_ind == k2p_ind
                            Extra_part_T2(i) = Extra_part_T2(i) + 1;
                            count_rho = count_rho+1;
                            coords_to_rho(count_rho,:) = [i,index_finder_rho(k1_ind,sz_rho_1+1)];
                            values_to_rho(count_rho) = -1;

                            count_rho = count_rho+1;
                            coords_to_rho(count_rho,:) = [i,index_finder_rho(k2_ind,sz_rho_2+1)];
                            values_to_rho(count_rho) = -1;
                        end
                    end

                    if k2_ind ~= k3p_ind && k2p_ind ~= k3_ind
                        Qi_T1 = SkkBZ(k1_ind,SkkBZ(k2_ind,k3p_ind));
                        parity = 1;
                        parity_p = 1;

                        if k2_ind > k3p_ind
                            temp = k2_ind;
                            k2_ind = k3p_ind;
                            k3p_ind = temp;
                            parity = -1;
                        end

                        if k2p_ind > k3_ind
                            temp = k2p_ind;
                            k2p_ind = k3_ind;
                            k3_ind = temp;
                            parity_p = -1;
                        end

                        count_T1 = count_T1 + 1;
                        coords_to_T1(count_T1,:) = [i,index_finder_T1(k1_k2_ordered_finder(k2_ind,k3p_ind),k1_k2_ordered_finder(k2p_ind,k3_ind),Qi_T1,Sz_T1)];
                        values_to_T1(count_T1) = -parity*parity_p;

                    end

                end

                if Sz_T2 == 3
                    Sz_T1 = 3;
                    if k2_ind > k1_ind && k1p_ind ~= k3_ind

                        Qi_T1 = SkkBZ(k1_ind,SkkBZ(k2_ind,k3p_ind));
                        parity_p = 1;

                        if k1p_ind > k3_ind
                            temp = k1p_ind;
                            k1p_ind = k3_ind;
                            k3_ind = temp;
                            parity_p = -1;
                        end

                        count_T1 = count_T1 + 1;
                        coords_to_T1(count_T1,:) = [i,index_finder_T1(k1_k2_ordered_finder(k1_ind,k2_ind),k1_k2_ordered_finder(k1p_ind,k3_ind),Qi_T1,Sz_T1)];
                        values_to_T1(count_T1) = +parity_p;
                    end
                end

                if Sz_T2 == 4
                    Sz_T1 = 3;
                    if k2p_ind > k1p_ind && k1_ind ~= k3p_ind

                        Qi_T1 = SkkBZ(k1_ind,SkkBZ(k2_ind,k3p_ind));
                        parity = 1;

                        if k1_ind > k3p_ind
                            temp = k1_ind;
                            k1_ind = k3p_ind;
                            k3p_ind = temp;
                            parity = -1;
                        end

                        count_T1 = count_T1 + 1;
                        coords_to_T1(count_T1,:) = [i,index_finder_T1(k1_k2_ordered_finder(k1_ind,k3p_ind),k1_k2_ordered_finder(k1p_ind,k2p_ind),Qi_T1,Sz_T1)];
                        values_to_T1(count_T1) = +parity;
                    end
                end

                if Sz_T2 == 7
                    Sz_T1 = 4;
                    if k2_ind > k1_ind && k1p_ind ~= k3_ind

                        Qi_T1 = SkkBZ(k1_ind,SkkBZ(k2_ind,k3p_ind));
                        parity_p = 1;

                        if k1p_ind > k3_ind
                            temp = k1p_ind;
                            k1p_ind = k3_ind;
                            k3_ind = temp;
                            parity_p = -1;
                        end

                        count_T1 = count_T1 + 1;
                        coords_to_T1(count_T1,:) = [i,index_finder_T1(k1_k2_ordered_finder(k1_ind,k2_ind),k1_k2_ordered_finder(k1p_ind,k3_ind),Qi_T1,Sz_T1)];
                        values_to_T1(count_T1) = +parity_p;
                    end
                end

                if Sz_T2 == 8
                    Sz_T1 = 4;
                    if k2p_ind > k1p_ind && k1_ind ~= k3p_ind

                        Qi_T1 = SkkBZ(k1_ind,SkkBZ(k2_ind,k3p_ind));
                        parity = 1;

                        if k1_ind > k3p_ind
                            temp = k1_ind;
                            k1_ind = k3p_ind;
                            k3p_ind = temp;
                            parity = -1;
                        end

                        count_T1 = count_T1 + 1;
                        coords_to_T1(count_T1,:) = [i,index_finder_T1(k1_k2_ordered_finder(k1_ind,k3p_ind),k1_k2_ordered_finder(k1p_ind,k2p_ind),Qi_T1,Sz_T1)];
                        values_to_T1(count_T1) = +parity;
                    end
                end

            end
        end
    end
end

coords_to_T1 = coords_to_T1(abs(values_to_T1)>0,:);
values_to_T1 = values_to_T1(abs(values_to_T1)>0);

M_T2_to_T1 = sparse(coords_to_T1(:,1),coords_to_T1(:,2),values_to_T1,dim_T2,dim_T1);

coords_to_D = coords_to_D(abs(values_to_D)>0,:);
values_to_D = values_to_D(abs(values_to_D)>0);

M_T2_to_D = sparse(coords_to_D(:,1),coords_to_D(:,2),values_to_D,dim_T2,dim_d);

coords_to_rho = coords_to_rho(abs(values_to_rho)>0,:);
values_to_rho = values_to_rho(abs(values_to_rho)>0);

M_T2_to_rho = sparse(coords_to_rho(:,1),coords_to_rho(:,2),values_to_rho,dim_T2,dim_rho);


[M_T1_to_D,M_T1_to_rho,Extra_part_T1] = T1_to_D(Nk,SkkBZ,DkkBZ);

M_T2_to_D = M_T2_to_D + M_T2_to_T1*M_T1_to_D;
M_T2_to_rho = M_T2_to_rho + M_T2_to_T1*M_T1_to_rho;
Extra_part_T2 = Extra_part_T2 + M_T2_to_T1*Extra_part_T1;

end

function [ksz_index,index_finder] = rho_index(Nk)
    num_sz = 2;
    ksz_index = zeros(num_sz*Nk,2); 
    index_finder = zeros(Nk,num_sz);
    count = 0;
    for sz = 0:num_sz-1
        for i = 1:Nk
            count = count + 1;
            ksz_index(count,:) = [i,sz];
            index_finder(i,sz+1) = count;
        end
    end
end


function [index,k1_k2_ordered,index_finder] = T1_index(Nk)
    num_Sz = 4;
    index = zeros(Nk*(Nk*(Nk-1)/2)^2*num_Sz,4);
    index_finder = zeros(Nk*(Nk-1)/2,Nk*(Nk-1)/2,Nk,num_Sz);
    k1_k2_ordered = zeros(Nk*(Nk-1)/2,2);

    count = 0;
    for k1 = 1:Nk-1
        for k2 = k1+1:Nk
            count = count + 1;
            k1_k2_ordered(count,:) = [k1,k2];
        end
    end
    assert(count==Nk*(Nk-1)/2)

    count = 0;
    for Sz = 1:num_Sz
        for Qi = 1:Nk
            for kp = 1:Nk*(Nk-1)/2
                for k = 1:Nk*(Nk-1)/2
                    count = count + 1;
                    index(count,:) = [k,kp,Qi,Sz];
                    index_finder(k,kp,Qi,Sz) = count;
                end
            end
        end
    end
end


function [M_T1_to_D,M_T1_to_rho,Extra_part] = T1_to_D(Nk,SkkBZ,DkkBZ)
    [index_d,index_finder_d] = D_indexing(Nk);
    [index_rho,index_finder_rho] = rho_index(Nk);
    [index_T1,T1_k1_k2_ordered,~] = T1_index(Nk);


    dim_T1 = length(index_T1);
    dim_d = length(index_d);
    dim_rho = length(index_rho);

    coords_to_D = zeros(dim_T1,2);
    values_to_D = zeros(dim_T1,1);
    coords_to_rho = zeros(dim_T1,2);
    values_to_rho = zeros(dim_T1,1);

    count_D = 0;
    count_rho = 0;


    Extra_part = zeros(dim_T1,1);

    for i = 1:dim_T1
        k = index_T1(i,1);
        kp = index_T1(i,2);
        Qi_T1 = index_T1(i,3);
        Sz_T1 = index_T1(i,4);
        k1_ind = T1_k1_k2_ordered(k,1);
        k2_ind = T1_k1_k2_ordered(k,2);
        k3_ind = DkkBZ(Qi_T1,SkkBZ(k1_ind,k2_ind));
        k1p_ind = T1_k1_k2_ordered(kp,1);
        k2p_ind = T1_k1_k2_ordered(kp,2);
        k3p_ind = DkkBZ(Qi_T1,SkkBZ(k1p_ind,k2p_ind));

        if Sz_T1 == 1
            Sz_d = 1;
            sz_rho = 0; 
                circlic_ind = [k1_ind,k2_ind,k3_ind;k2_ind,k3_ind,k1_ind;k3_ind,k1_ind,k2_ind];
                circlic_p_ind = [k1p_ind,k2p_ind,k3p_ind;k2p_ind,k3p_ind,k1p_ind;k3p_ind,k1p_ind,k2p_ind];

                for or = 1:3
                    for or_p = 1:3
                        k1_ind = circlic_ind(or,1);
                        k2_ind = circlic_ind(or,2);
                        k3_ind = circlic_ind(or,3);
                        k1p_ind = circlic_p_ind(or_p,1);
                        k2p_ind = circlic_p_ind(or_p,2);
                        k3p_ind = circlic_p_ind(or_p,3);

                        if (k1_ind==k1p_ind && k2_ind==k2p_ind && k3_ind==k3p_ind) || ...
                                (k1_ind==k1p_ind && k2_ind==k3p_ind && k3_ind==k2p_ind)

                            temp = (k1_ind==k1p_ind)*(k2_ind==k2p_ind)*(k3_ind==k3p_ind)...
                                -(k1_ind==k1p_ind)*(k2_ind==k3p_ind)*(k3_ind==k2p_ind);


                            Extra_part(i) = Extra_part(i) + 1/3*temp;

                            count_rho = count_rho+1;
                            coords_to_rho(count_rho,:) = [i,index_finder_rho(k1_ind,sz_rho+1)];
                            values_to_rho(count_rho) = -1/3*temp;

                            count_rho = count_rho+1;
                            coords_to_rho(count_rho,:) = [i,index_finder_rho(k2_ind,sz_rho+1)];
                            values_to_rho(count_rho) = -1/3*temp;

                            count_rho = count_rho+1;
                            coords_to_rho(count_rho,:) = [i,index_finder_rho(k3_ind,sz_rho+1)];
                            values_to_rho(count_rho) = -1/3*temp;

                        end

                        if k1_ind == k1p_ind
                            Qi_d = SkkBZ(k2_ind,k3_ind);
                            assert(SkkBZ(k2p_ind,k3p_ind)==Qi_d)
                            parity = 1;
                            parity_p = 1;

                            if k2_ind > k3_ind
                                temp = k2_ind;
                                k2_ind = k3_ind;
                                k3_ind = temp;
                                parity = -1;
                            end

                            if k2p_ind > k3p_ind
                                temp = k2p_ind;
                                k2p_ind = k3p_ind;
                                k3p_ind = temp;
                                parity_p = -1;
                            end

                            count_D = count_D + 1;
                            coords_to_D(count_D,:) = [i,index_finder_d(k2_ind,k2p_ind,Qi_d,Sz_d)];
                            values_to_D(count_D) = parity_p*parity;

                        end
                    end
                end
        end


        if Sz_T1 == 2
            Sz_d = 2;
            sz_rho = 1;
                circlic_ind = [k1_ind,k2_ind,k3_ind;k2_ind,k3_ind,k1_ind;k3_ind,k1_ind,k2_ind];
                circlic_p_ind = [k1p_ind,k2p_ind,k3p_ind;k2p_ind,k3p_ind,k1p_ind;k3p_ind,k1p_ind,k2p_ind];

                for or = 1:3
                    for or_p = 1:3
                        k1_ind = circlic_ind(or,1);
                        k2_ind = circlic_ind(or,2);
                        k3_ind = circlic_ind(or,3);
                        k1p_ind = circlic_p_ind(or_p,1);
                        k2p_ind = circlic_p_ind(or_p,2);
                        k3p_ind = circlic_p_ind(or_p,3);

                        if (k1_ind==k1p_ind && k2_ind==k2p_ind && k3_ind==k3p_ind) || ...
                                (k1_ind==k1p_ind && k2_ind==k3p_ind && k3_ind==k2p_ind)

                            temp = 1/3*((k1_ind==k1p_ind)*(k2_ind==k2p_ind)*(k3_ind==k3p_ind)...
                                -(k1_ind==k1p_ind)*(k2_ind==k3p_ind)*(k3_ind==k2p_ind));

                            Extra_part(i) = Extra_part(i) + temp;

                            count_rho = count_rho+1;
                            coords_to_rho(count_rho,:) = [i,index_finder_rho(k1_ind,sz_rho+1)];
                            values_to_rho(count_rho) = -temp;

                            count_rho = count_rho+1;
                            coords_to_rho(count_rho,:) = [i,index_finder_rho(k2_ind,sz_rho+1)];
                            values_to_rho(count_rho) = -temp;

                            count_rho = count_rho+1;
                            coords_to_rho(count_rho,:) = [i,index_finder_rho(k3_ind,sz_rho+1)];
                            values_to_rho(count_rho) = -temp;

                        end

                        if k1_ind == k1p_ind
                            Qi_d = SkkBZ(k2_ind,k3_ind);
                            assert(SkkBZ(k2p_ind,k3p_ind)==Qi_d)
                            parity = 1;
                            parity_p = 1;

                            if k2_ind > k3_ind
                                temp = k2_ind;
                                k2_ind = k3_ind;
                                k3_ind = temp;
                                parity = -1;
                            end

                            if k2p_ind > k3p_ind
                                temp = k2p_ind;
                                k2p_ind = k3p_ind;
                                k3p_ind = temp;
                                parity_p = -1;
                            end

                            count_D = count_D + 1;
                            coords_to_D(count_D,:) = [i,index_finder_d(k2_ind,k2p_ind,Qi_d,Sz_d)];
                            values_to_D(count_D) = parity_p*parity;

                        end
                    end
                end
        end

        if Sz_T1 == 3
            if k3_ind == k3p_ind
                Sz_d = 1;
                assert(k1_ind<k2_ind && k1p_ind<k2p_ind)
                Qi_d = SkkBZ(k1_ind,k2_ind);
                assert(SkkBZ(k1p_ind,k2p_ind)==Qi_d)
                count_D = count_D + 1;
                coords_to_D(count_D,:) = [i,index_finder_d(k1_ind,k1p_ind,Qi_d,Sz_d)];
                values_to_D(count_D) = 1;
            end

            circlic_12 = [k1_ind,k2_ind;k2_ind,k1_ind];
            circlic_p_12 = [k1p_ind,k2p_ind;k2p_ind,k1p_ind];
            sign = [1,-1];

            for or = 1:2
                for or_p = 1:2
                    sign_or = sign(or);
                    sign_or_p = sign(or_p);
                    k1_ind = circlic_12(or,1);
                    k2_ind = circlic_12(or,2);
                    k1p_ind = circlic_p_12(or_p,1);
                    k2p_ind = circlic_p_12(or_p,2);

                    if k1_ind == k1p_ind
                        Sz_d = 3;
                        Qi_d = SkkBZ(k2_ind,k3_ind);
                        count_D = count_D + 1;
                        coords_to_D(count_D,:) = [i,index_finder_d(k2_ind,k2p_ind,Qi_d,Sz_d)];
                        values_to_D(count_D) = sign_or*sign_or_p;
                    end

                    if k1_ind==k1p_ind && k2_ind==k2p_ind && k3_ind==k3p_ind
                        temp = 1/2;

                        Extra_part(i) = Extra_part(i) + temp*sign_or*sign_or_p;

                        count_rho = count_rho+1;
                        coords_to_rho(count_rho,:) = [i,index_finder_rho(k1_ind,1)];
                        values_to_rho(count_rho) = -temp*sign_or*sign_or_p;

                        count_rho = count_rho+1;
                        coords_to_rho(count_rho,:) = [i,index_finder_rho(k2_ind,1)];
                        values_to_rho(count_rho) = -temp*sign_or*sign_or_p;

                        count_rho = count_rho+1;
                        coords_to_rho(count_rho,:) = [i,index_finder_rho(k3_ind,2)];
                        values_to_rho(count_rho) = -temp*sign_or*sign_or_p;


                    end

                end
            end


        end

        if Sz_T1 == 4
            if k3_ind == k3p_ind
                Sz_d = 2;
                assert(k1_ind<k2_ind && k1p_ind<k2p_ind)
                Qi_d = SkkBZ(k1_ind,k2_ind);
                assert(SkkBZ(k1p_ind,k2p_ind)==Qi_d)
                count_D = count_D + 1;
                coords_to_D(count_D,:) = [i,index_finder_d(k1_ind,k1p_ind,Qi_d,Sz_d)];
                values_to_D(count_D) = 1;
            end

            circlic_12 = [k1_ind,k2_ind;k2_ind,k1_ind];
            circlic_p_12 = [k1p_ind,k2p_ind;k2p_ind,k1p_ind];
            sign = [1,-1];

            for or = 1:2
                for or_p = 1:2
                    sign_or = sign(or);
                    sign_or_p = sign(or_p);
                    k1_ind = circlic_12(or,1);
                    k2_ind = circlic_12(or,2);
                    k1p_ind = circlic_p_12(or_p,1);
                    k2p_ind = circlic_p_12(or_p,2);

                    if k1_ind == k1p_ind
                        Sz_d = 3;
                        Qi_d = SkkBZ(k2_ind,k3_ind);
                        count_D = count_D + 1;
                        coords_to_D(count_D,:) = [i,index_finder_d(k3_ind,k3p_ind,Qi_d,Sz_d)];
                        values_to_D(count_D) = sign_or*sign_or_p;
                    end

                    if (k1_ind==k1p_ind && k2_ind==k2p_ind && k3_ind==k3p_ind)
                        temp = 1/2;

                        Extra_part(i) = Extra_part(i) + temp*sign_or*sign_or_p;

                        count_rho = count_rho+1;
                        coords_to_rho(count_rho,:) = [i,index_finder_rho(k1_ind,2)];
                        values_to_rho(count_rho) = -temp*sign_or*sign_or_p;

                        count_rho = count_rho+1;
                        coords_to_rho(count_rho,:) = [i,index_finder_rho(k2_ind,2)];
                        values_to_rho(count_rho) = -temp*sign_or*sign_or_p;

                        count_rho = count_rho+1;
                        coords_to_rho(count_rho,:) = [i,index_finder_rho(k3_ind,1)];
                        values_to_rho(count_rho) = -temp*sign_or*sign_or_p;


                    end

                end
            end


        end
    end

    coords_to_D = coords_to_D(abs(values_to_D)>0,:);
    values_to_D = values_to_D(abs(values_to_D)>0);

    M_T1_to_D = sparse(coords_to_D(:,1),coords_to_D(:,2),values_to_D,dim_T1,dim_d);

    coords_to_rho = coords_to_rho(abs(values_to_rho)>0,:);
    values_to_rho = values_to_rho(abs(values_to_rho)>0);

    M_T1_to_rho = sparse(coords_to_rho(:,1),coords_to_rho(:,2),values_to_rho,dim_T1,dim_rho);

end
