function [x_pred, V, G, P, th]=rkiteration(A,B,C,D,V,tau,c,x,y)
%RKIteration robust Kalman iteration 
%
%   [x_pred, V_next, G]=RKIteration(A,B,C,D,V,tau,c,x,y) performs one 
%   iteration of the "delayed" robust Kalman estimator for the nominal 
%    discrete-time model
%      
%      x[n+1] = Ax[n] + Bw[n]         {State equation}
%        y[n] = Cx[n] + Dw[n]         {Measurements}
%
%   with disturbance w with variance I.
%   The equation of the robust "delayed" estimator:
%
%      x[n+1|n] = Ax[n|n-1] + G (y[n] - Cx[n|n-1])
%
%   where y[n] past measurement, x[n|n-1] past estimate, and G robust
%   Kalman gain. 
%
%   RKALMAN returns the new estimate x[n+1|n], the least feavorable
%   covariance matrix V_next of the estimation error x[n+1]-x[n+1|n], and 
%   the estimator gain G.
%
%   The robust estimator is designed knowing that the true model belongs to 
%   a ball about the nominal one. The models in that ball are such that 
%   the Tau-divergence between them and the nominal model is less than a 
%   certain tolerance. To design the estimator gain G, it is necessary to 
%   specify: 
%    
%    - the tolerance c (striclty positive)
%    - the parameter tau (in the interval [0,1]) of the Tau-divergence
%    - the least favorable covariance matrix V of the previous estimation  
%      error x[n]-x[n|n-1].
%
%   For more details see: "Robust Kalman filtering under incremental 
%   model perturbations" by M. Zorzi
%
%   See also RKALMAN, MAXTOL.


%   Author(s): Mattia Zorzi 20-8-2015


% parameters 

n=size(A,1);
p=size(C,1);
m=size(B,2);

Q=B*B';
R=D*D';

% robust Kalman gain
G = A*V*C'*(C*V*C'+R)^-1;

% compute the the prediction
x_pred=(A*x'+G*(y-C*x'))';


% update the least favorable covariance matrix V
P = (A-G*C)*V*(A-G*C)'+(B-G*D)*(B-G*D)';
L=chol(P)';
value=1;
t1=0;
if tau==1
    t2=10/max(eig(P)); % because exp(700) is (more or less) the maximum number that matlab can store
else
    e = eig(P);
    r = max(abs(e));
    t2=(1-10^-5)*((1-tau)*r)^-1;
end
while abs(value)>=10^-9
   th=0.5*(t1+t2);
   if tau==0
       value=trace(inv(eye(n)-th*P)-eye(n)) + log(det(eye(n)-th*P))-c;
   end
   if tau>0 & tau<1
      value=trace(-1/(tau*(1-tau))*(eye(n)-(1-tau)*th*L'*L)^(tau/(tau-1))+1/(1-tau)*(eye(n)-(1-tau)*th*L'*L)^(1/(tau-1))+1/tau*eye(n))-c;
   end
   if tau==1
       value=trace(th*L'*L*expm(th*L'*L)-expm(th*L'*L)+eye(n))-c;
   end
   if value>0
        t2=th;
   else
        t1=th;
   end
end
Vold=V;
if tau==1
    V = L*expm(th*L'*L)*L';
else
    V = L*(eye(n)-(1-tau)*th*L'*L)^(1/(tau-1))*L';
end


