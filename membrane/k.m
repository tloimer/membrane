function k=k(T,a)
%K(T,A)     Effective thermal conductivity of the mixture-filled membrane [W/mK].
%
%  Calls EPSILON, KG, KL.

% epsilon: Lückengrad - volume fraction of pores
% epsilon=0.6;
km=23; % ceramic
%km = 0.22; % Celgard

k = (1-epsilon)*km + epsilon*(a.*kg(T)+(1-a).*kl(T));
