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