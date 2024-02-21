function [theta,gamma,Yhat,sse,lamc,gcv,bic]=TPST_est(B,Q2,K,Y,lambda)

n = length(Y);
W = B*Q2;
WW = W'*W;
WY = W'*Y;
D = Q2'*K*Q2;
nq = size(Q2,2);
Iq = eye(nq);
[~,flag] = chol(WW);
if(flag~=0)
   WW = WW+diag(ones(nq,1)*1e-12);
end
Ainv = chol(WW,'upper');
A = Ainv'\Iq;
ADA = A*D*A';
C = eig(full(ADA));
    
nl = length(lambda);
sse_all = []; df_all = []; gcv_all = []; bic_all = []; 
theta_all = []; gamma_all = [];
for il=1:nl
    Lam = lambda(il);
    Dlam = Lam*D;
    Winv = (WW+Dlam)\Iq;
%     if(flag~=0)
%         Hmtx=W*Winv*W';
%     end
    theta = Winv*WY;
    theta_all = [theta_all,theta];
    gamma = Q2*theta;
    gamma_all = [gamma_all,gamma];
    Yhat = W*theta;
    sse = sum((Y-Yhat).^2);
    sse_all = [sse_all,sse];
    df = sum(1./(1+C*Lam));
%     if(flag==0)
%         df=sum(1./(1+C*Lam));
%     end
%     if(flag~=0)
%         df=trace(Hmtx);
%     end
    df_all = [df_all,df];
    gcv = n*sse/(n-df)^2;
    gcv_all = [gcv_all,gcv];   
    bic = log(sse/n)+df*log(n)/n;
    bic_all = [bic_all,bic];
end
[gcv,lam_ind] = min(gcv_all);
lamc = lambda(lam_ind);
theta = theta_all(:,lam_ind);
gamma = gamma_all(:,lam_ind);
df = df_all(lam_ind);
sse = sse_all(lam_ind);
bic = bic_all(lam_ind);
Yhat = W*theta;
