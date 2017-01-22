function [x, G, V, P, th]=rkalman(A,B,C,D,y,c,tau)
%RKALMAN  robust Kalman estimator.
%
%   [x_e,G,V] = RKALMAN(A,B,C,D,y,c,tau) compute the "delayed" robust 
%   Kalman estimator for the nominal discrete-time model
%      
%      x[n+1] = Ax[n] + Bv[n]         {State equation}
%        y[n] = Cx[n] + Dv[n]         {Measurements}
%
%   with disturbance and measurement noise v with variance I. The robust 
%   "delayed" estimator uses only past measurements up to y[n-1] to 
%   generate the optimal robust estimate x_e[n] of x[n] and is easier to 
%   embed in digital control loops. The equation of the robust "delayed" 
%   estimator:
%
%      x_e[n+1|n] = Ax_e[n|n-1] + G (y[n] - Cx_e[n|n-1])
%
%   RKALMAN returns the delayed estimate x_e, the estimator gain G, the 
%   least feavorable covariance matrix V of the estimation error 
%   x[n+1]-x_e[n+1|n].
%
%   The robust estimator is designed knowing that the true model belongs to 
%   a ball about the nominal one. The models in that ball are such that 
%   the Tau-divergence between them and the nominal model is less than a 
%   certain tolerance. To design the robust Kalman estimator, it is 
%   necessary to specify: 
%    
%    - the tolerance c (striclty positive)
%    - the parameter tau (in the interval [0,1]) of the Tau-divergence
%    - the measurements y
%
%   For more details see: "Robust Kalman filtering under incremental 
%   model perturbations" by M. Zorzi
%
%   See also RKITERATION, MAXTOL.


%   Author(s): Mattia Zorzi 20-8-2015



% check the inputs
if nargin==6
    tau=1;
end


% parameters 
n=size(A,1);
p=size(C,1);
m=size(B,2);
T=size(y,1);

% transform the model
A=A-B*D'*(D*D')^-1*C;
B=[(B*(eye(m)-D'*(D*D')^-1*D)*B')^0.5 zeros(n,m-n)];
D=[zeros(p,n) (D*D')^0.5];
Q=B*B';
R=D*D';

% check conditions
assert(c>0,'Tolerance c must be positive')
assert(rank(ctrb(A,Q))>=n,'The model must be reachable')
assert(rank(ctrb(A',C'))>=n,'The model must be observable')
cN = maxtol(A,B,C,D,tau,2*n);
if c>cN
    warning('Tolerance c is too large: the filter gain may not exist')
end

% init
x=zeros(T,n);
V=zeros(n,n,T+1);
P=zeros(n,n,T+1);
G=zeros(n,p,T+1);
th=zeros(T,1);
V(:,:,1)=eye(n);


% iterative part
for k=1:T
    [x(k+1,:), V(:,:,k+1), G(:,:,k+1), P(:,:,k+1), th(k)]=rkiteration(A,B,C,D,V(:,:,k),tau,c,x(k,:),y(k,:));
end

% resize
x=x(2:T+1,:);
P=P(:,:,2:T+1);
V=V(:,:,2:T+1);
G=G(:,:,2:T+1);


