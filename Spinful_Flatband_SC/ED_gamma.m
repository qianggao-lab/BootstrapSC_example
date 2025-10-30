function [stiffness_ED,runtime] = ED_gamma(Nx,Ny,N,A)
checker(Nx,Ny,N,A)
format long
kx0 = -1/2+1/Nx/2:1/Nx:1/2;
ky0 = -1/2+1/Ny/2:1/Ny:1/2;

a1 = [1, 0];
a2 = [0, 1];

A_mat = [a1; a2];
B_mat = 2*pi*inv(A_mat)';
b1 = B_mat(1,:);
b2 = B_mat(2,:);


b_mat = [b1;b2];

[KX, KY] = meshgrid(kx0, ky0);
k0_ind   = [KX(:), KY(:)];
k0 = k0_ind*b_mat;

Nk = length(k0_ind);

SkkBZ = zeros(Nk,Nk);
for i = 1:Nk
    for j = 1:Nk
        temp1 = mod(k0_ind(i,1)+k0_ind(j,1),1);
        if temp1 >= 1/2 
            temp1 = temp1 - 1; 
        end
        temp2 = mod(k0_ind(i,2)+k0_ind(j,2),1);
        if temp2 >= 1/2 
            temp2 = temp2 - 1; 
        end
        SkkBZ(i,j) = Ny*(temp1*Nx+(Nx+1)/2-1) + (temp2*Ny+(Ny+1)/2);
    end
end
SkkBZ = round(SkkBZ);
DkkBZ = flip(SkkBZ,2);

F_Q_kkp_ssp_A = Hamiltonian(Nk,A,k0,DkkBZ);
assert(mod(N,2)==0)
tic
[~, E0, ~] = flatband_ED_gamma(F_Q_kkp_ssp_A, SkkBZ, DkkBZ, N, N/2);
runtime = toc;
stiffness_ED = E0/((A(1))^2+A(2)^2)/(Nx*Ny)^2/4/4;
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

function [F_Q_kkp_ssp_A] = Hamiltonian(Nk,A,k0,DkkBZ)
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


F_Q_kkp_ssp_A = zeros(Nk,Nk,Nk,2,2);


for k1 = 1:Nk
    for k1p = 1:Nk
        for Qi = 1:Nk
            for sz = 1:2
                for szp = 1:2

                    xi_sz = (-1)^(sz-1);
                    xi_szp = (-1)^(szp-1);

                    k2 = DkkBZ(k1,Qi);
                    k2p = DkkBZ(k1p,Qi);

                    for alpha = 1:2
                        for beta = 1:2
                            for pm = 1:2 
                                for I = 1:4
                                    F_Q_kkp_ssp_A(k1,k1p,Qi,sz,szp) = F_Q_kkp_ssp_A(k1,k1p,Qi,sz,szp) + 1/2*xi_sz*xi_szp...
                                        *D_AB_k(DkkBZ(k1,k2),I,pm,alpha)*D_AB_k(DkkBZ(k2p,k1p),I,pm,beta)...
                                        *(conj(U_k_alpha_flat_spin(alpha,k1,sz))*U_k_alpha_flat_spin(alpha,k2,sz)...
                                        *conj(U_k_alpha_flat_spin(beta,k2p,szp))*U_k_alpha_flat_spin(beta,k1p,szp));
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
end

function [H, E0, info] = flatband_ED_gamma(F, Skk, Dkk, N, Nup)
Nk = size(Skk,1);
L  = 2*Nk;
assert(ndims(F)==5 && all([size(F,1),size(F,2),size(F,3)] == Nk) ...
    && size(F,4)==2 && size(F,5)==2, 'F must be [Nk,Nk,Nk,2,2].');
assert(isequal(size(Skk),[Nk,Nk]) && isequal(size(Dkk),[Nk,Nk]), 'Skk/Dkk size mismatch.');
assert(L <= 64, 'This fast path requires 2*Nk <= 64 (L=%d).', L);

tol     = 1e-10;
maxit   = 1000;

Ndown = N - Nup;
assert(Nup>=0 && Ndown>=0 && Nup<=Nk && Ndown<=Nk, 'Invalid (Nup,Ndown).');

Gamma = int32(Dkk(1,1));

Q_list = 1:Nk;
Q_list = int32(Q_list);

orb_k = int32([(1:Nk), (1:Nk)]);
orb_s = int32([ones(1,Nk), 2*ones(1,Nk)]);
ks2orb = @(k,s) int32((s==1)*k + (s==2)*(Nk+k));

orb_plusQ  = zeros(L, Nk, 'int32');
orb_minusQ = zeros(L, Nk, 'int32');
for o = 1:L
    k = orb_k(o); s = orb_s(o);
    for Q = 1:Nk
        kp = Skk(k, Q);
        km = Dkk(k, Q);
        orb_plusQ(o,  Q) = ks2orb(kp, s);
        orb_minusQ(o, Q) = ks2orb(km, s);
    end
end


preMask = zeros(1,L,'uint64');
for pos = 1:L
    if pos==1, preMask(pos) = uint64(0);
    else,      preMask(pos) = bitshift(uint64(1), pos-1) - 1; end
end

F11 = F(:,:,:,1,1);
F12 = F(:,:,:,1,2);
F21 = F(:,:,:,2,1);
F22 = F(:,:,:,2,2);

[UpByK]   = spin_basis_byK(Nk, Nup,  Skk, Gamma);
[DownByK] = spin_basis_byK(Nk, Ndown, Skk, Gamma);
[basis64, occCache] = combine_spin_bases_to_Gamma(UpByK, DownByK, Dkk, Gamma, Nk, Nup, Ndown);
nStates = numel(basis64);
if nStates == 0
    warning('Empty Γ sector for the chosen (Nup,Ndown).'); E0 = NaN;
    info = struct('Nk',Nk,'N',N,'Nup',Nup,'Ndown',Ndown,'nStates',0,'nnzH',0,'Gamma',double(Gamma),'converged',false);
    return;
end

[basis_sorted, ord] = sort(basis64);
find_idx = @(st) ord(binary_search_uint64(basis_sorted, st));

est_per_row = max(1, N*N*double(numel(Q_list)));
nnz_est = max(1, int64(ceil(0.5 * est_per_row * nStates)));


I = zeros(nnz_est,1,'uint32'); J = zeros(nnz_est,1,'uint32');
V = complex(zeros(nnz_est,1), zeros(nnz_est,1));
p = uint32(0);
for idx = uint32(1):uint32(nStates)
    st  = basis64(idx);
    occ = occCache(idx,:);
    for t4 = 1:N
        o4 = occ(t4);  s4 = orb_s(o4);  k4 = orb_k(o4);
        [st1, sgn4, ok4] = apply_annihil64(st, o4, preMask);
        if ~ok4, continue; end
        for q = 1:numel(Q_list)
            Q = Q_list(q);
            o3 = orb_minusQ(o4, Q);
            if bitget(st1, o3), continue; end
            [st2, sgn3] = apply_create64(st1, o3, preMask);
            if o3==o4
                occ2 = occ;
            else
                occ2 = occ; occ2(t4) = o3; occ2 = sort(occ2);
            end
            for t2 = 1:N
                o2 = occ2(t2); s2 = orb_s(o2);
                [st3, sgn2, ok2] = apply_annihil64(st2, o2, preMask);
                if ~ok2, continue; end
                o1 = orb_plusQ(o2, Q);
                if bitget(st3, o1), continue; end
                [st4, sgn1] = apply_create64(st3, o1, preMask);
                jdx = find_idx(st4);
                if jdx < idx, continue; end

                k1 = orb_k(o1);
                if s2==1
                    if s4==1, g = F11(k1,k4,Q); else, g = F12(k1,k4,Q); end
                else
                    if s4==1, g = F21(k1,k4,Q); else, g = F22(k1,k4,Q); end
                end
                amp = (sgn4 * sgn3 * sgn2 * sgn1) * g;

                p = p + 1;
                if p > numel(I)
                    grow = max(1, ceil(0.5*double(p)));
                    I = [I; zeros(grow,1,'uint32')];
                    J = [J; zeros(grow,1,'uint32')];
                    V = [V; complex(zeros(grow,1),zeros(grow,1))];
                end
                I(p) = idx; J(p) = jdx; V(p) = V(p) + amp;
            end
        end
    end
end
I = I(1:p); J = J(1:p); V = V(1:p);

H = sparse(double(I), double(J), V, double(nStates), double(nStates));

H = H + H' - spdiags(diag(H), 0, nStates, nStates);



nnzH = nnz(H);
eigs_opts = struct('tol', tol, 'maxit', maxit);

[~, D, flag] = eigs(-H, 1, 'largestreal', eigs_opts);
E0 = -real(D(1,1));
converged = (flag == 0);

info = struct('Nk',Nk,'N',N,'Nup',Nup,'Ndown',Ndown,'nStates',nStates,'nnzH',nnzH,'Gamma',double(Gamma),'converged',converged);
end

function C = spin_basis_byK(Nk, Nspin, Skk, Kzero)
C = cell(1, Nk);
for K=1:Nk, C{K} = struct('masks', uint64([]), 'occs', zeros(0, Nspin, 'int32')); end
if Nspin == 0
    C{Kzero}.masks = uint64(0); C{Kzero}.occs = zeros(1,0,'int32'); return;
end
comb = int32(1:Nspin); M = int32(Nk); done = false;
while ~done
    Ksum = Kzero;
    for t=1:Nspin, Ksum = Skk(Ksum, comb(t)); end
    mask = uint64(0);
    for t=1:Nspin, mask = bitset(mask, comb(t), 1); end
    bucket = C{Ksum};
    bucket.masks(end+1,1) = mask;
    bucket.occs(end+1,:)   = comb(:).';
    C{Ksum} = bucket;
    [comb, done] = next_combination(comb, M);
end
end

function [basis64, occCache] = combine_spin_bases_to_Gamma(UpByK, DownByK, Dkk, Gamma, Nk, Nup, Ndown)
total = 0;
for Ku = 1:Nk
    Kd = Dkk(Gamma, Ku);
    total = total + numel(UpByK{Ku}.masks) * numel(DownByK{Kd}.masks);
end
basis64  = zeros(total,1,'uint64');
occCache = zeros(total, Nup+Ndown, 'int32');
p = 0;
for Ku = 1:Nk
    Kd   = Dkk(Gamma, Ku);
    U    = UpByK{Ku};  D = DownByK{Kd};
    nu   = numel(U.masks);  nd = numel(D.masks);
    if nu==0 || nd==0, continue; end
    for iu = 1:nu
        up_mask = U.masks(iu);
        up_occ  = U.occs(iu,:);
        for id = 1:nd
            dn_mask = D.masks(id);
            dn_occ  = D.occs(id,:);
            st = bitor(up_mask, bitshift(dn_mask, Nk));
            p = p + 1;
            basis64(p) = st;
            occCache(p,:) = int32([up_occ, Nk + dn_occ]);
        end
    end
end
end

function [c, done] = next_combination(c, M)
N = numel(c); done = false;
i = N;
while i >= 1 && c(i) == M - (N - i), i = i - 1; end
if i == 0, done = true; return; end
c(i) = c(i) + 1;
for j = i+1:N, c(j) = c(j-1) + 1; end
end

function pos = binary_search_uint64(arr, key)
lo = int32(1); hi = int32(numel(arr));
while lo <= hi
    mid = bitshift(lo + hi, -1);
    v = arr(mid);
    if v < key
        lo = mid + 1;
    elseif v > key
        hi = mid - 1;
    else
        pos = mid; return;
    end
end
error('State not found in basis (internal bug).');
end

function [st_out, sgn, ok] = apply_annihil64(st_in, pos, preMask)
if ~bitget(st_in, pos), st_out = st_in; sgn = 0; ok = false; return; end
before = popcount64(bitand(st_in, preMask(pos)));
sgn = sign_from_parity(before);
st_out = bitset(st_in, pos, 0); ok = true;
end

function [st_out, sgn] = apply_create64(st_in, pos, preMask)
before = popcount64(bitand(st_in, preMask(pos)));
sgn = sign_from_parity(before);
st_out = bitset(st_in, pos, 1);
end

function s = sign_from_parity(n)
if bitand(n,1)==0, s = +1; else, s = -1; end
end

function n = popcount64(u)
persistent LUT
if isempty(LUT), LUT = uint8(sum(dec2bin(0:255)=='1',2)); end
b = typecast(u, 'uint8');
n = sum(double(LUT(double(b)+1)));
end