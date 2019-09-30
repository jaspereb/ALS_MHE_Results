function [P_est, Q_est, R_est]=MEALS(Y,Me,N,k0,AK,Abar,Ck,Gk,Hk)
% This is the Memory Efficient ALS code used for the journal paper
% experiments
disp("Should MEALS be enforcing diagonality?");

if ~iscell(Ck)
    tild_C=kron(speye(Me),Ck);
    [LP,LN]=size(Ck);
else
    tild_C=sparse(blkdiag(Ck{1,k0:k0+Me-1}));
    [LP,LN]=size(Ck{2});
end

if ~iscell(Gk)
    tild_G=kron(speye(Me),Gk);    [LN,LR]=size(Gk);
else
    tild_G=sparse(blkdiag(Gk{1,k0-1:k0+Me-2}));
    [LN,LR]=size(Gk{2});
end
[~,LQ] = size(Hk);

tild_A=kron(speye(Me),speye(LN))-[[sparse(LN,(Me-1)*LN);sparse(blkdiag(Abar{1,k0:k0+Me-2}))] sparse(LN*Me,LN)];
tild_AK=sparse(blkdiag(AK{1,k0-1:k0+Me-2}));  %Need to add H here if not eye

CA_inv=tild_C/tild_A; %tild V
tild_B=CA_inv*tild_G;
tild_D=CA_inv*tild_AK;   %Negative appears when calculating AK in 'Pre-calculate for ALS' section
tild_M=aux_mat(LP,LP*N,1);


% tild_H=[speye(Ly), sparse(Ly,(N-1)*Ly)];

A_cal1_cell=cell(Me-N+1,1);  A_cal2_cell=cell(Me-N+1,1);  A_cal3_cell=cell(Me-N+1,1);
D_n=duplic_mat(LN);   D_q=duplic_mat(LQ);   D_r=duplic_mat(LR);

for i=0:Me-N %Can use parfor
    tild_E=[repmat(Abar{k0-1},1,(i+1));sparse((Me-1)*LN,(i+1)*LN)]; %Added (i+1) to column dim, had to widen Abar
    tild_F=CA_inv*tild_E;
    
    %Auxiliary selection matrices
    tild_Si=aux_mat(LP*N,LP*Me,LP*i+1);
    tild_Ji=(aux_mat(LR*(i+1),LR*Me,1))';
    tild_Ui=(aux_mat(LQ*(i+1),LQ*Me,1))';
    tild_Pi=aux_mat(LP*(N-1),LP*Me,(LP*(i+1)+1));
    tild_Oi=(aux_mat(LQ,LQ*Me,(LQ*(i+1)+1)))';
    
    Gamma_i = tild_Si*tild_F;
    Gamma_bar = tild_M*Gamma_i;
    
    Omega_i=full(tild_Si*tild_B*tild_Ji);
    Omega_bar=full(tild_M*Omega_i);
    
    Phi_i=tild_Si*tild_D*tild_Ui;
    Phi_bar=full(tild_M*Phi_i);
    
    Psi_i=[Hk',(tild_Pi*tild_D*tild_Oi)']';
    
    
    %============================ Memory Efficient A_cal =======================================
    P_term = zeros(N*LP*LP, LN*LN);
    Q_term = zeros(N*LP*LP, LR*LR);
    R_term = zeros(N*LP*LP, LQ*LQ);
    
    for j = 1:(i+1)
        z_gamma = Decomp(LN,j,i);
        P_term = P_term + (kron((Gamma_bar*z_gamma),(Gamma_i*z_gamma)));
        
        z_omega = Decomp(LR,j,i);
        Q_term = Q_term + (kron((Omega_bar*z_omega),(Omega_i*z_omega)));
        
        z_phi = Decomp(LQ,j,i);
        R_term = R_term + (kron((Phi_bar*z_phi),(Phi_i*z_phi)));
    end
 
    A_cal1_cell{i+1,1}=P_term*D_n;
    A_cal2_cell{i+1,1}=Q_term*D_r;
    A_cal3_cell{i+1,1}=(R_term+kron(Hk,Psi_i))*D_q;
    
end
disp('Part 1 finished');

%============================== Calculate R_bar ==================================
Rbar = [];
idxesl = [];
idxesr = [];

for col = 1:(Me-N+1) %Is the i value
    for row = 1:N
        leftIdx = col + row + k0;
        rightIdx = 1 + col + k0;
        idxesl(row,col) = leftIdx;
        idxesr(row,col) = rightIdx;
        
        zz = Y{leftIdx}*(Y{rightIdx}');
        zz = zz(:);
        
        Rbar = [Rbar; zz];
    end
end
b = Rbar;

%============================== SDP code ==================================
A_cal=cell2mat([A_cal1_cell, A_cal2_cell, A_cal3_cell]);

cvx_begin sdp
cvx_solver Mosek

variable P(LN,LN) symmetric semidefinite diagonal
variable Q(LR,LR) symmetric semidefinite diagonal
variable R(LP,LP) symmetric semidefinite diagonal

Pdum = P';
mask = triu(true(LN,LN),0)';
P0bar = Pdum(mask);

Qdum = Q';
mask = triu(true(LR,LR),0)';
Qbar = Qdum(mask);

Rdum = R';
mask = triu(true(LQ,LQ),0)';
Rbar = Rdum(mask);

X=[P0bar;Qbar;Rbar];

minimize( sum((A_cal*X-b).^2));
subject to
    Q(1,1)==Q(3,3);
    Q(2,2)==Q(4,4);
    Q(7,7)==Q(9,9);
    Q(8,8)==Q(10,10);
    Q(5,5)==Q(3,3);
    Q(2,2)==Q(6,6);
    Q(7,7)==Q(11,11);
    Q(8,8)==Q(12,12);
    
cvx_end

P_est = P;
Q_est = Q;
R_est = R;
disp('Solver Completed');

end

function PM=permut_mat(p,N)
a=[speye(p);sparse(p*N-p,p)];
a=kron(speye(p),a);
b=[a;sparse(p,p^2)];
b=kron(ones(N,1),b);
PM=b(1:(p*N)^2,:);
end

function DM=duplic_mat(N)
a=1:(N*(N+1)/2);
b=tril(ones(N));
b(b==1)=a;
b=b+triu(b',1);
index=b(:);
DM=sparse(1:N^2,index,ones(N^2,1));
end

function AM=aux_mat(r,c,l)
%Make an auxiliary matrix (script M)
AM=[sparse(r,(l-1)),speye(r),sparse(r,(c-r-l+1))];
end
