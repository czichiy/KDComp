function KD=CalculateKD(Dic, Nd, sigma, Ruu)

% Matrix form calculation  (fast)
KD = zeros(Nd,Nd,Nd^2);
cst =  1/sqrt(det(eye(size(Ruu))+4/sigma^2*Ruu));
M = inv(eye(size(Ruu))+sigma^2/4*inv(Ruu));
MDic = M * Dic;
DicMDic = Dic'*MDic;
DiagDicMDic = sum(Dic.*MDic);
P2lp =  ones(Nd,1)*DiagDicMDic + DiagDicMDic'*ones(1,Nd)+ 2*DicMDic;


for i = 1 : Nd
    for j = 1 : Nd
          P1 = ones(Nd,1)*sum(Dic.^2) + sum(Dic.^2)'*ones(1,Nd) + Dic(:,i)'*Dic(:,i) + Dic(:,j)'*Dic(:,j);
          P1 = -4*P1;
          P2ij = DicMDic(i,i) + DicMDic(j,j)+2*DicMDic(i,j);
          uij = Dic(:,i) + Dic(:,j);
          P2lpij =  2*(MDic'* (uij*ones(1,Nd)) + (ones(Nd,1)*uij')*MDic);
          P2 = P2ij + P2lp + P2lpij;
          P = (P1 + P2)/(8*sigma^2);
          KD(:,:,(i-1)*Nd+j)=cst*exp(P);
    end
end



% 
 %% equivalent component-wise form (easy to read)
% for i = 1 : Nd
%     for j = 1 : Nd
%         for l = 1 : Nd
%             for p = 1 : Nd
%                 ui = Dic(:,i);
%                 uj = Dic(:,j);
%                 ul = Dic(:,l);
%                 up = Dic(:,p);
%                 KDij= exp(-(4*(ui'*ui+uj'*uj+up'*up+ul'*ul) + (ui+uj+ul+up)'*M*(ui+uj+ul+up))/(8*sigma^2));
%             end
%         end
%         KD(:,:,(i-1)*Nd+j) = cst*KDij;
%     end
% end
