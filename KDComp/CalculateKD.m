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
        for l = 1 : Nd
            for p = 1 : Nd
                ui = Dic(:,i);
                uj = Dic(:,j);
                ul = Dic(:,l);
                up = Dic(:,p);
                KDij= exp(-(4*(ui'*ui+uj'*uj+up'*up+ul'*ul) + (ui+uj+ul+up)'*M*(ui+uj+ul+up))/(8*sigma^2));
            end
        end
        KD(:,:,(i-1)*Nd+j) = cst*KDij;
    end
end
