function RkkD=CalculateRkkD(Dic, Nd, sigma, Ruu)

% Matrix form calculation  (fast)
P1 = ones(Nd,1)*sum(Dic.^2) + sum(Dic.^2)'*ones(1,Nd);
P1 = -2*P1;
M = inv(eye(size(Ruu))+sigma^2/2*inv(Ruu));
MDic = M * Dic;
P2 = ones(Nd,1)*sum(Dic.*MDic) + sum(Dic.*MDic)'*ones(1,Nd)+ 2*Dic'*MDic;
Rkk1 = (P1+P2)/(4*sigma^2);
RkkD = 1/sqrt(det(eye(size(Ruu))+2/sigma^2*Ruu))*exp(Rkk1);