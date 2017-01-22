function cN=maxtol(A,B,C,D,tau,N)
%MAXTOL maximum tolerance for which the robust Kalman estimator is 
%       guaranteed to converge.
%
%   [c] = MAXTOL(A,B,C,D,tau) returns the maximum tolerance for which the 
%   robust Kalman estimator is guaranteed to converge. 
%
%   The inputs of MAXTOL are the nominal discrete-time model associated to
%   the robust Kalman filter
%      
%      x[n+1] = Ax[n] + Bv[n]         {State equation}
%        y[n] = Cx[n] + Dv[n]         {Measurements}
%
%   with disturbance and measurement noise v with variance I (independent), 
%   the parameter tau (in [0,1]) of the Tau-divergence.
%
%   For more details see: "Robust Kalman filtering under incremental 
%   model perturbations" by M. Zorzi 
%
%   See also RKALMAN, RKITERATION.

%   Author(s): Mattia Zorzi 20-8-2015


% check inputs
if nargin==5
    N=2*n;
end

n=size(A,1);
m=size(B,2);
p=size(C,1);

% construction of the Gramians
DR=kron(eye(N),D);
Re=[];
Ob=[];
ObR=[];
H=[];
L=[];
for k=1:N
    Re=[Re A^(k-1)*B]; 
    Ob=[C*A^(k-1); Ob];
    ObR=[A^(k-1); ObR];
    T=[];
    for l=1:N
        if l<=N-(k-1)
            T=[T zeros(p,m)];
        else
            T=[T C*A^(l-N+(k-1)-1)*B];
        end
    end
    H=[T; H];
    T=[];
    for l=1:N
        if l<=N-(k-1)
            T=[T zeros(n,m)];
        else
            T=[T A^(l-N+(k-1)-1)*B];
        end
    end
    L=[T; L];
end

Om = Ob'*inv(DR*DR'+ H*H')*Ob;
J = ObR -L*H'*inv(DR*DR'+H*H')*Ob;
M= L*inv(eye(N*m)+H'*inv(DR*DR')*H)*L';

% computation of phiN tilde
phiNtilde = 1/max(eig(M));

% computation of phiN 
value=1;
t1=0;
t2=(1-10^-10)*phiNtilde;
while abs(value)>=10^-9
       theta=0.5*(t1+t2);
       Sth = -eye(N*n)+theta*M;
       Omth = Om +J'*theta*(Sth)^-1*J;
       value=-min(eig(Omth));
       if value>0
            t2=theta;
       else
            t1=theta;
       end
end
phiN=min(phiNtilde,theta);

% computation of the initial condition (that is P_q bar) 
Pq=dare(A',C',B*B',D*D');
lambda=min(eig(Pq));

% computation of theta_N
if tau==1
    thN=-log(1-phiN*lambda)/lambda;
else
    thN=(1-(1-lambda*phiN)^(1-tau))/((1-tau)*lambda);
end
if tau>=0 & tau<1
    thNmax=((1-tau)*max(eig(Pq)))^-1;
else
    thNmax=+inf;
end
% computation of cN

if thN>=thNmax
    cN=+inf;
else
    if tau==0
        cN=(log(det(eye(n)-thN*Pq))+trace((eye(n)-thN*Pq)^-1)-n);
    end
    if tau>0 & tau<1
        Lq=chol(Pq)';
        cN=trace(-1/(tau*(1-tau))*(eye(n)-thN*(1-tau)*Lq'*Lq)^(tau/(tau-1))+1/(1-tau)*(eye(n)-thN*(1-tau)*Lq'*Lq)^(1/(tau-1))+1/tau*eye(n));
    end
    if tau==1
        Lq=chol(Pq)';
        cN=trace(expm(thN*Lq'*Lq)*(thN*Lq'*Lq-eye(n))+eye(n));
    end
end




