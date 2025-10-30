%{
Copyright (c) 2025 Qiang Gao
SPDX-License-Identifier: MIT
Author: Qiang Gao <gq201277@gmail.com>
Description: Bootstrap the stiffness of a chiral flatband superconductor
%}

function [nu,stiffness] = stiffness_chiralSC(Nx,Ny,M,s,plotResults)
checker(Nx,Ny,M,s);
N_list = 2:2:Nx*Ny;
stiffness_x = zeros(size(N_list)); 
stiffness_y = zeros(size(N_list)); 

for i = 1:length(N_list)
    fprintf('Running... Step %d/%d\t',N_list(i)/2,Nx*Ny/2)
    tic
    stiffness_x(i) = bootstrap_chiralSC(Nx,Ny,N_list(i),[0.002,0],M,s);
    stiffness_y(i) = bootstrap_chiralSC(Nx,Ny,N_list(i),[0,0.002],M,s);
    toc
end

[nu,stiffness] = gceconvertion(Nx,Ny,N_list,M,s,stiffness_x,stiffness_y);

if plotResults
    plot_results(nu,stiffness,M,s)
end
end

function checker(Nx,Ny,M,s)
if ~(isfinite(Nx) && (Nx == fix(Nx))) || ~(isfinite(Ny) && (Ny == fix(Ny)))  || Nx<=0 || Ny <=0
    error('Invalid inputs: Nx and Ny must be positive integers')
elseif Nx~=Ny
    error('This algorithm only works for Nx=Ny')
elseif mod(Nx,2)~=0
    error('This algorithm only works for even by even lattices')
end

if ~(isfinite(M) && (M == fix(M))) || M<=0
    error('Invalid inputs: M must be a positive integer')
elseif M > 1
    warning('For M>1, you should choose a much larger lattice.')
end

if s <= 0
    error('Invalid inputs: s must be a positive real number')
elseif s < 1 
    warning('For s<1, you should choose a larger lattice.')
end
end

function [nu,stiffness] = gceconvertion(Nx,Ny,N_list,M,s,stiffness_x,stiffness_y)
if ~isequal(sort(N_list),(2:2:Nx*Ny))
    error('The list of particle numbers must include all even particle sectors from 2 to Nx*Ny')
end
kx0 = -1/2+1/Nx:1/Nx:1/2+1/Nx/2;
ky0 = -1/2+1/Ny:1/Ny:1/2+1/Ny/2;

b1 = 2*pi*[1,1/sqrt(3)];
b2 = 2*pi*[1,-1/sqrt(3)];

B_mat = [b1;b2];

[KX, KY] = meshgrid(kx0, ky0);
k0_ind   = [KX(:), KY(:)];
k0 = k0_ind*B_mat;

Nk = length(k0_ind);

omega1 = (b1(1) + 1i*b1(2))/2;
omega2 = (b2(1) + 1i*b2(2))/2;

k0 = k0 +[-2*pi/Nx,0];

k0_complex = k0(:,1)+1i*k0(:,2);


periodic_zeta = zeta_haldane(k0_complex, omega2, omega1);

normalization = sqrt(1 + abs(s*periodic_zeta).^(2*M));

u_A = (s*periodic_zeta).^M./normalization;
u_B = 1./normalization;

f_k = sort(abs(u_A./u_B).^2);

f_k_symm = f_k(1:2:end-1);

e_f_k = elemSymAll(f_k_symm);

N_list = [0,N_list];
stiffness_x = [0,stiffness_x];
stiffness_y = [0,stiffness_y];


alpha = -3:0.005:5;

nu_alpha_xx = zeros(size(alpha));
Ds_alpha_xx = zeros(size(alpha));

for ai = 1:length(alpha)
    nu_alpha_xx(ai) = sum(exp(-2*alpha(ai)*N_list).*N_list.*e_f_k)/sum(exp(-2*alpha(ai)*N_list).*e_f_k)/Nk;
    Ds_alpha_xx(ai) = sum(exp(-2*alpha(ai)*N_list).*stiffness_x.*e_f_k)/sum(exp(-2*alpha(ai)*N_list).*e_f_k);
end

nu_alpha_yy = zeros(size(alpha));
Ds_alpha_yy = zeros(size(alpha));

for ai = 1:length(alpha)
    nu_alpha_yy(ai) = sum(exp(-2*alpha(ai)*N_list).*N_list.*e_f_k)/sum(exp(-2*alpha(ai)*N_list).*e_f_k)/Nk;
    Ds_alpha_yy(ai) = sum(exp(-2*alpha(ai)*N_list).*stiffness_y.*e_f_k)/sum(exp(-2*alpha(ai)*N_list).*e_f_k);
end

nu_alpha = (nu_alpha_xx+nu_alpha_yy)/2;
Ds_alpha = sqrt(Ds_alpha_xx.*Ds_alpha_yy);


nu_selected = 1:-0.05:0;
idx = interp1(nu_alpha, 1:numel(nu_alpha),nu_selected, 'nearest', 'extrap');
nu = nu_alpha(idx);
stiffness = Ds_alpha(idx);
end

function e = elemSymAll(x)          
n  = numel(x);
e  = zeros(1,n+1);  e(1) = 1;      
for xi = x(:).'                   
    e(2:end) = e(2:end) + xi*e(1:end-1);  
end
end

function f = zeta_haldane(Z, omega1, omega2)
L1     = 2*omega1;
L2     = 2*omega2;
A_cell = abs(imag(conj(L1).*L2)); 
B      = -pi ./ A_cell;    

wweier = 'zeta';
f = arrayfun(@(z) ellipWeier(wweier,omega1, omega2, z).value+conj(z)*B, Z);
end


function stiffness = bootstrap_chiralSC(Nx,Ny,N,A0,M,s)
format long
assert(Nx==Ny)
assert(mod(Nx,2)==0)
assert(mod(N,2)==0)

kx0 = -1/2+1/Nx:1/Nx:1/2+1/Nx/2;
ky0 = -1/2+1/Ny:1/Ny:1/2+1/Ny/2;


b1 = 2*pi*[1,1/sqrt(3)];
b2 = 2*pi*[1,-1/sqrt(3)];

B_mat = [b1;b2];

[KX, KY] = meshgrid(kx0, ky0);
k0_ind   = [KX(:), KY(:)];

k0 = k0_ind*B_mat;

Nk = length(k0_ind);

[M_ex,mask,M_cont,M_T_to_2,M_T_to_1,Extra_part,Dkk,M_k_Qmk] = con_mat(k0_ind,N);

if norm(A0)<0.00001
    error("The inserted flat gauge field A cannot be too small")
end

A = A0+[-2*pi/Nx,0];
omega1 = (b1(1) + 1i*b1(2))/2;
omega2 = (b2(1) + 1i*b2(2))/2;
[F_kkp,F_k] = Hamiltonian_chiral_SC(Nk,A,k0,Dkk,M,s,omega1,omega2);


ops = sdpsettings('solver','mosek','verbose',0);

M_D = sdpvar(Nk,Nk,Nk,'hermitian','complex');
rho = sdpvar(Nk,1,'full','real');

M_D = mask.*M_D;
Cons = [rho == reshape(M_cont*reshape(M_D,[],1),[],1),sum(rho) == N];

M_G = reshape(M_ex*reshape(M_D,[],1),Nk,Nk,Nk);

M_T = reshape(M_T_to_2*reshape(M_D,[],1)+M_T_to_1*reshape(rho,[],1)+Extra_part,[Nk*(Nk-1)/2,Nk*(Nk-1)/2,Nk]);

M_Q = permute(M_D,[2,1,3]);
for Qi = 1:Nk
    M_Q(:,:,Qi) = M_Q(:,:,Qi)...
        +mask(:,:,Qi).*diag(ones(Nk,1)-rho-M_k_Qmk(:,:,Qi)*rho);
end


for Qi = 1:Nk
    Cons = [Cons, rho(Qi) >= 0, M_D(:,:,Qi)>=0, (M_G(:,:,Qi)+M_G(:,:,Qi)')/2 + diag(rho)>=0,...
        M_Q(:,:,Qi)>=0,M_T(:,:,Qi)+M_T(:,:,Qi)'>=0];
end

Objective = 100*real(transpose(F_kkp)*reshape(M_G(:,:,:),[],1)+ reshape(F_k,1,[])*reshape(rho,[],1));

diagnose = optimize(Cons, Objective, ops);

M_g_value = value(M_G);
rho_value = value(rho);
Energy = transpose(F_kkp)*reshape(M_g_value(:,:,:),[],1)+ reshape(F_k,1,[])*reshape(rho_value,[],1);


yalmip('clear')
factor = 2/sqrt(3);
if diagnose.problem ~= 0 
    fprintf('The solver did not return a converged result. Please be careful about this data point: Nx=%d, Ny=%d, N=%d',Nx,Ny,N)
end

stiffness = factor*Energy/Nk/Nk/(A0(1)^2+A0(2)^2)/2;
assert(max(abs(imag(stiffness)))<1e-10);
stiffness = real(stiffness);
end

function [F_kkp,F_k] = Hamiltonian_chiral_SC(Nk,A,k0,DkkBZ,M,s,omega1,omega2)
k0_pA = k0 + A;
k0_pA_complexified = k0_pA(:,1)+1i*k0_pA(:,2);

periodic_zeta = zeta_haldane(k0_pA_complexified, omega2, omega1);

normalization = sqrt(1 + abs(s*periodic_zeta).^(2*M));

u_A_pA = (s*periodic_zeta).^M./normalization;
u_B_pA = 1./normalization;

k_index_ordered_F = zeros(Nk^3,3);

count = 0;
for Qi = 1:Nk
    for kp = 1:Nk
        for k = 1:Nk
            count = count + 1;
            k_index_ordered_F(count,:) = [k,kp,Qi];
        end
    end
end

F_kkp = zeros(length(k_index_ordered_F),1);
F_k = zeros(Nk,1);
for i = 1:length(F_kkp)
    k1 = k_index_ordered_F(i,1);
    k1p = k_index_ordered_F(i,2);
    Qi = k_index_ordered_F(i,3);
    k2 = DkkBZ(k1,Qi);
    k2p = DkkBZ(k1p,Qi);
    F_kkp(i) = conj(u_B_pA(k1))*u_A_pA(k2)*conj(u_A_pA(k2p))*u_B_pA(k1p);
end

for i = 1:Nk
    k1 = i;
    k1p = i;
    temp = 0;
    for Qi = 1:Nk
        k2 = DkkBZ(k1,Qi);
        k2p = k2;
        temp = temp + conj(u_B_pA(k1))*u_A_pA(k2)*conj(u_A_pA(k2p))*u_B_pA(k1p);
    end
    F_k(i) = temp;
end
end

function [M_ex,mask,M_tr,M_T2D,M_T2rho,Extra,Dkk,M_k_Qmk] = con_mat(k0_ind,N)
Nk = length(k0_ind);

[Skk, Dkk, mask, M_k_Qmk] = lattice(k0_ind,Nk);

M_ex = G2D(Nk,Skk,Dkk);

M_tr = D2rho(Nk,N,Dkk);

[M_T2D,M_T2rho,Extra] = Con_T2(Nk,Skk,Dkk);

end

function [Skk, Dkk, mask, M_k_Qmk] = lattice(k0_ind,Nk)
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

mask = ones(Nk,Nk,Nk);
for Qi = 1:Nk
    for k = 1:Nk
        for kp = 1:Nk
            if Dkk(Qi,k) <= k || Dkk(Qi,kp) <= kp
                mask(k,kp,Qi) = 0;
            end
        end
    end
end

M_k_Qmk = zeros(Nk,Nk,Nk);
for Qi = 1:Nk
    for k = 1:Nk
        M_k_Qmk(k,Dkk(Qi,k),Qi) = 1;
    end
end
end

function M_ex = G2D(Nk,SkkBZ,DkkBZ)
    [index_d,index_finder_d] = D_indexing(Nk);
    index_g = G_indexing(Nk);

    dim_d = length(index_d);
    dim_g = length(index_g);

    coords = zeros(dim_g,2);
    values = zeros(dim_g,1);

    count = 0;

    for i = 1:dim_g
        k1 = index_g(i,1);
        k1p = index_g(i,2);
        Qi = index_g(i,3);

        k2 = DkkBZ(k1,Qi);
        k2p = DkkBZ(k1p,Qi);

        k1t = k1;
        k2t = k2p;
        k1tp = k1p;
        k2tp = k2;

        assert(SkkBZ(k2t,k1t) == SkkBZ(k2tp,k1tp))
        Qi_d = SkkBZ(k2t,k1t);
            if k2t ~= k1t && k2tp ~= k1tp

                es = 1;
                es_p = 1;

                if k1t > k2t
                    k1t = k2t;
                    es = -1;
                end

                if k1tp > k2tp
                    k1tp = k2tp;
                    es_p = -1;
                end

                count = count + 1;
                coords(count,:) = [i,index_finder_d(k1t,k1tp,Qi_d)];
                values(count) = -es*es_p;
            end

    end

    coords = coords(abs(values)>0,:);
    values = values(abs(values)>0);

    M_ex = sparse(coords(:,1),coords(:,2),values,dim_g,dim_d);
end

function M_cont = D2rho(Nk,N,DkkBZ)
[index_d,index_finder_d] = D_indexing(Nk);
dim_d = length(index_d);

dim_rho = Nk;

coords = zeros(6*Nk*dim_rho,2);
values = zeros(6*Nk*dim_rho,1);

count = 0;

for k = 1:dim_rho
    for Qt = 1:Nk
        k1t = k;
        k2t = DkkBZ(Qt,k);
        k1tp = k;
        k2tp = DkkBZ(Qt,k);

        if k1t ~= k2t && k1tp ~= k2tp

            es = 1;
            es_p = 1;

            if k1t > k2t
                k1t = k2t;
                es = -1;
            end

            if k1tp > k2tp
                k1tp = k2tp;
                es_p = -1;
            end

            count = count + 1;
            coords(count,:) = [k,index_finder_d(k1t,k1tp,Qt)];
            values(count) = es*es_p/(N-1);
        end

    end

end

coords = coords(abs(values)>0,:);
values = values(abs(values)>0);
M_cont = sparse(coords(:,1),coords(:,2),values,dim_rho,dim_d);

end

function [M_T2_to_D,M_T2_to_rho,Extra_part_T2] = Con_T2(Nk,SkkBZ,DkkBZ)
[index_d,index_finder_d] = D_indexing(Nk);
[index_rho,index_finder_rho] = rho_index(Nk);
[index_T1,T1_k1_k2_ordered,index_finder_T1] = T1_index(Nk);
[index_T2,T2_k1_k2,~] = T2_index(Nk);

dim_T1 = length(index_T1);
dim_d = length(index_d);
dim_rho = length(index_rho);

coords_to_D = zeros(dim_T1,2);
values_to_D = zeros(dim_T1,1);
coords_to_rho = zeros(dim_T1,2);
values_to_rho = zeros(dim_T1,1);

count_D = 0;
count_rho = 0;


Extra_part_T1 = zeros(dim_T1,1);

for i = 1:dim_T1
    k = index_T1(i,1);
    kp = index_T1(i,2);
    Qi_T1 = index_T1(i,3);
    k1_ind = T1_k1_k2_ordered(k,1);
    k2_ind = T1_k1_k2_ordered(k,2);
    k3_ind = DkkBZ(Qi_T1,SkkBZ(k1_ind,k2_ind));
    k1p_ind = T1_k1_k2_ordered(kp,1);
    k2p_ind = T1_k1_k2_ordered(kp,2);
    k3p_ind = DkkBZ(Qi_T1,SkkBZ(k1p_ind,k2p_ind));


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

                Extra_part_T1(i) = Extra_part_T1(i) + 1/3*temp;

                count_rho = count_rho+1;
                coords_to_rho(count_rho,:) = [i,index_finder_rho(k1_ind)];
                values_to_rho(count_rho) = -1/3*temp;

                count_rho = count_rho+1;
                coords_to_rho(count_rho,:) = [i,index_finder_rho(k2_ind)];
                values_to_rho(count_rho) = -1/3*temp;

                count_rho = count_rho+1;
                coords_to_rho(count_rho,:) = [i,index_finder_rho(k3_ind)];
                values_to_rho(count_rho) = -1/3*temp;

            end

            if k1_ind == k1p_ind
                Qi_d = SkkBZ(k2_ind,k3_ind);
                assert(SkkBZ(k2p_ind,k3p_ind)==Qi_d)
                es = 1;
                es_p = 1;

                if k2_ind > k3_ind
                    k2_ind = k3_ind;
                    es = -1;
                end

                if k2p_ind > k3p_ind
                    k2p_ind = k3p_ind;
                    es_p = -1;
                end

                count_D = count_D + 1;
                coords_to_D(count_D,:) = [i,index_finder_d(k2_ind,k2p_ind,Qi_d)];
                values_to_D(count_D) = es_p*es;

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

Extra_part_T1 = sparse(Extra_part_T1);


k1_k2_ordered_finder = zeros(Nk,Nk);

count = 0;
for k1 = 1:Nk-1
    for k2 = k1+1:Nk
        count = count + 1;
        k1_k2_ordered_finder(k1,k2) = count;
    end
end
assert(count==Nk*(Nk-1)/2)

dim_T2 = length(index_T2);
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

for i = 1:dim_T2
    k = index_T2(i,1);
    kp = index_T2(i,2);
    Qi_T2 = index_T2(i,3);

    k1_ind = T2_k1_k2(k,1);
    k2_ind = T2_k1_k2(k,2);
    k3_ind = DkkBZ(SkkBZ(k1_ind,k2_ind),Qi_T2);
    k1p_ind = T2_k1_k2(kp,1);
    k2p_ind = T2_k1_k2(kp,2);
    k3p_ind = DkkBZ(SkkBZ(k1p_ind,k2p_ind),Qi_T2);

    assert( k2_ind > k1_ind && k2p_ind > k1p_ind)
        if k3_ind == k3p_ind
            Qi_d = SkkBZ(k1_ind,k2_ind);

            count_D = count_D + 1;
            coords_to_D(count_D,:) = [i,index_finder_d(k1_ind,k1p_ind,Qi_d)];
            values_to_D(count_D) = 2;

            if k1_ind == k1p_ind && k2_ind == k2p_ind
                Extra_part_T2(i) = Extra_part_T2(i) + 1;
                count_rho = count_rho+1;
                coords_to_rho(count_rho,:) = [i,index_finder_rho(k1_ind)];
                values_to_rho(count_rho) = -1;

                count_rho = count_rho+1;
                coords_to_rho(count_rho,:) = [i,index_finder_rho(k2_ind)];
                values_to_rho(count_rho) = -1;
            end

        end

        if k3p_ind ~= k1_ind && k3p_ind ~= k2_ind && k3_ind ~= k1p_ind && k3_ind ~= k2p_ind
            Qi_T1 = SkkBZ(k1_ind,SkkBZ(k2_ind,k3p_ind));
            es = 1;
            es_p = 1;

            if k3p_ind < k1_ind
                k2_ind = k1_ind;
                k1_ind = k3p_ind;
                es = 1;
            elseif k3p_ind < k2_ind
                k2_ind = k3p_ind;
                es = -1;
            end

            if k3_ind < k1p_ind
                k2p_ind = k1p_ind;
                k1p_ind = k3_ind;
                es_p = 1;
            elseif k3_ind < k2p_ind
                k2p_ind = k3_ind;
                es_p = -1;
            end

            count_T1 = count_T1 + 1;
            coords_to_T1(count_T1,:) = [i,index_finder_T1(k1_k2_ordered_finder(k1_ind,k2_ind),k1_k2_ordered_finder(k1p_ind,k2p_ind),Qi_T1)];
            values_to_T1(count_T1) = -es*es_p;
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

Extra_part_T2 = sparse(Extra_part_T2);

M_T2_to_D = M_T2_to_D + M_T2_to_T1*M_T1_to_D;
M_T2_to_rho = M_T2_to_rho + M_T2_to_T1*M_T1_to_rho;
Extra_part_T2 = Extra_part_T2 + M_T2_to_T1*Extra_part_T1;
end

function [ksz_index,index_finder] = rho_index(Nk)
    ksz_index = zeros(Nk,1); 
    index_finder = zeros(Nk,1);
    count = 0;
        for i = 1:Nk
            count = count + 1;
            ksz_index(count,:) = i;
            index_finder(i) = count;
        end
end

function [index,index_finder] = D_indexing(Nk)
k_index_ordered = zeros(Nk^3,3);
index_finder = zeros(Nk,Nk,Nk);

count = 0;
for Qi = 1:Nk
    for kp = 1:Nk
        for k = 1:Nk
            count = count + 1;
            k_index_ordered(count,:) = [k,kp,Qi];
            index_finder(k,kp,Qi) = count;
        end
    end
end
index = k_index_ordered;
end

function [index_g,index_finder_g] = G_indexing(Nk)
    k_index_ordered = zeros(Nk^3,3);
    index_finder = zeros(Nk,Nk,Nk);

    count = 0;
        for Qi = 1:Nk
            for kp = 1:Nk
                for k = 1:Nk
                    count = count + 1;
                    k_index_ordered(count,:) = [k,kp,Qi];
                    index_finder(k,kp,Qi) = count;
                end
            end
        end

    index_g = k_index_ordered;

    if nargout>=2
        index_finder_g = index_finder;
    end
end

function [index,k1_k2_ordered,index_finder] = T1_index(Nk)
index = zeros(Nk*(Nk*(Nk-1)/2)^2,3);
index_finder = zeros(Nk*(Nk-1)/2,Nk*(Nk-1)/2,Nk);
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
for Qi = 1:Nk
    for kp = 1:Nk*(Nk-1)/2
        for k = 1:Nk*(Nk-1)/2
            count = count + 1;
            index(count,:) = [k,kp,Qi];
            index_finder(k,kp,Qi) = count;
        end
    end
end
end

function [index,k1_k2,index_finder] = T2_index(Nk)
index = zeros(Nk*(Nk*(Nk-1)/2)^2,3);
index_finder = zeros(Nk*(Nk-1)/2,Nk*(Nk-1)/2,Nk);
k1_k2 = zeros(Nk*(Nk-1)/2,2);

count = 0;
for k1 = 1:Nk-1
    for k2 = k1+1:Nk
        count = count + 1;
        k1_k2(count,:) = [k1,k2];
    end
end

count = 0;
for Qi = 1:Nk
    for kp = 1:Nk*(Nk-1)/2
        for k = 1:Nk*(Nk-1)/2
            count = count + 1;
            index(count,:) = [k,kp,Qi];
            index_finder(k,kp,Qi) = count;
        end
    end
end
end

function plot_results(nu,stiffness_GCE,M,s)
figure('Color','w');
hold on


colororder({'#0072BD','#D95319'})

plot(nu, stiffness_GCE, '-s', ...
    'LineWidth',1.6, ...
    'MarkerSize',8, ...
    'MarkerFaceColor','w', ...
    'DisplayName','Bootstrap');


ax = gca;
ax.LineWidth   = 1.2;
ax.FontName    = 'Times New Roman';
ax.FontSize    = 14;
ax.TickDir     = 'out';
ax.Box         = 'on';
ax.GridAlpha   = 0.3; 
ax.GridLineStyle = '--';
grid on

xlabel('$\nu$', ...
    'Interpreter','latex', ...
    'FontSize',16, ...
    'FontName','Times New Roman');

ylabel('$D_s/V$', ...
    'Interpreter','latex', ...
    'FontSize',16, ...
    'FontName','Times New Roman');

tstr = sprintf('$A=%d, B=%d, M=%d, s=%d$',1, 2, M, s);
title(tstr, ...
    'Interpreter','latex', ...
    'FontSize',18, ...
    'FontWeight','bold', ...
    'FontName','Times New Roman');

leg = legend('Location','northeast');
set(leg, ...
    'Interpreter','latex', ...
    'FontSize',12, ...
    'Box','off');

hold off
end