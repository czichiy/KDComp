%{
Code Flow Logic

1. Initialization
   ├── Read input parameters: RLZ, eta, sigma, std_x, std_n, N, rho, a_lin, Nd
   ├── Create dictionary: Dic
   ├── Initialize variables: e2, e, Ed2, Ekd, RR

2. Perform RLZ filter operations
   ├── Loop: rlz = 1 to RLZ
   │   ├── Generate input signal: x
   │   ├── Calculate output signal: y
   │   ├── Add noise: dn
   │   ├── Initialize weights: alpha
   │   ├── Loop: n = 1 to N
   │   │   ├── Compute kernel function: kx
   │   │   ├── Estimate error: d_est
   │   │   ├── Update error: e(n)
   │   │   ├── Update weights: alpha
   │   └── Update error sum of squares: e2

3. Theoretical model computation
   ├── Preliminary computation: Ruu, RkkD, KD, KDvec
   ├── Optimal solution computation: alpha_optD, JminD
   ├── Initialization: CvD
   ├── Transient error computation
   │   ├── Loop: n = 1 to N
   │   │   ├── Compute TD
   │   │   ├── Update CvD
   │   │   ├── Compute theoretical error: e2_th

4. Compute steady-state error: JEMSED_infty, JMSED_infty

5. Plot results
   ├── Plot simulation curve: sim_curv
   ├── Plot theoretical error curves: e2_th, JMSED_infty


%}


% KLMS model using Gaussian kernel with uniform fixed dictionary.
% 1) IA input / temporal correlated input 
% 2) MC simulations and the model for MSE
%
% running time: apprx. 1min with 100 realization and 1E4 points.
%
% Contact: JIN CHIY 

tic
clc, clear

RLZ =100;         % number of realizations
eta = 0.05;       % step size
sigma=0.25;       % kernel bandwidth


std_x = 0.5;      % input std
std_n = 0.05;     % noise std

N = 30000;         % sequence length

rho = 0.5;        % input correaltion

a_lin = [0.5, -0.3]';    % linear transformation coefficients

Nd = 16;          % Dictionary size
Dic1 = [-1: 2/(sqrt(Nd)-1) :1];
Dic = [kron(ones(1,sqrt(Nd)),Dic1); kron(Dic1,ones(1,sqrt(Nd)))];     % uniform grid dictionary

% initializations
e2 = zeros(N,1);
e = zeros(N,1);
Ed2 = 0;
Ekd = zeros(Nd,1);
RR = 0;

%filtering
for rlz = 1 : RLZ
     if mod(rlz,20)==1, rlz, end

    % --input signal (IA assumptoin)  
%    x = randn(N,1);     
%    x = std_x*[x,rho*x+sqrt(1-rho^2)*randn(N,1)];
    % -- input signal (corr in time)
     x = std_x*filter(1,[1,-rho],sqrt(1-rho^2)*randn(N+50,1));
     x = x(50:end);
     x = [x(1:end-1),x(2:end)];
     
     %dic = std_x*filter(1,[1,-rho],sqrt(1-rho^2)*randn(Nd+50,1));
     %dic = x(50:50+Nd);
     %Dic = [dic(1:end-1);dic(2:end)];
     
    
    y = x*a_lin;               % linear output
    d =  y - 0.5*(y.^2)+0.1*y.^3;      % nonlinear output
    noise = std_n*randn(N,1);
    dn = d + noise;

    alpha = zeros(Nd,1);
    for n = 1 : N
        un = x(n,:)';
        kx = exp(-1/2/sigma^2*(un'*un+sum(Dic.^2) - 2*un'*Dic))';
        Ekd = Ekd+dn(n)*kx;
        d_est= alpha'*kx;
        e(n) = dn(n) - d_est;
        alpha = alpha + eta * e(n)*kx;
    end
    e2 = e2 + e.^2;
    Ed2 = Ed2 + mean(dn.^2);
end
e2 = e2/RLZ;
Ekd = Ekd/RLZ/N;
Ed2 = Ed2/RLZ;




%% Theoretical model 


%-- preliminary computation ...
%Ruu = std_x^2*eye(2);
Ruu = std_x^2*[1,rho;rho,1];

RkkD = CalculateRkkD(Dic, Nd, sigma, Ruu);
KD=CalculateKD(Dic, Nd, sigma, Ruu);
KDvec=reshape(KD,Nd^2,[]);

%-- Optimum
alpha_optD = inv(RkkD)*Ekd;
JminD = Ed2 - Ekd'*alpha_optD;

%-- Initialization
CvD = alpha_optD*alpha_optD';



%-- Transient error
for n = 1 : N
       if mod(n,1000) == 1, n, end
        TD = CalculateTD(KDvec,CvD, Nd);
        CvD = CvD - eta * (RkkD*CvD + CvD*RkkD) + eta^2*TD+eta^2*RkkD*JminD;
        e2_th(n) = JminD + trace(RkkD*CvD);
end
toc

%-- Steady state error
G = eye(Nd^2) - eta*(kron(eye(Nd),RkkD)+kron(RkkD,eye(Nd))) + eta^2*KDvec;
rkkD = reshape(RkkD,[],1);
JEMSED_infty=eta^2*JminD*rkkD'*inv(eye(Nd^2)-G)*rkkD;
JMSED_infty = JminD + JEMSED_infty;

sim_curv = filter(hamming(10)/sum(hamming(10)),1,e2);
figure, semilogy(sim_curv(10:end)), hold on
semilogy(e2), hold on
semilogy(e2_th, 'r', 'linewidth', 3), hold on
semilogy([1:N],JMSED_infty*ones(N,1),'m-.','linewidth',3)
legend('Simulated curve','Theoretical transient MSE', 'Theoretical steady-state MSE')
xlabel('Iteration $n$','fontsize',14,'interpreter', 'latex')
ylabel('$J_{\rm{MSE},\omega}$','fontsize',14,'interpreter', 'latex')
ylim([0.002,0.1])
grid on