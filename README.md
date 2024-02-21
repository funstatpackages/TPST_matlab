# TPST_matlab
This repo contains the implementation code for the article "Nonparametric Regression for 3D Point Cloud Learning."

1. Simulation Main Functions:
main_eg1_SP.m: main file for simulation study in Section 6.1 with domain Omega_1.
main_eg1_HS.m: main file for simulation study in Section 6.1 with domain Omega_2.


2. Data Generator:
dataGeneratorRandomSP.m: generating simulated data (random design) for simulation study in Section 6.1 with domain Omega_1.
dataGeneratorRandomHS.m: generating simulated data (random design) for simulation study in Section 6.1 with domain Omega_2.
popGeneratorSP.m: generating true data over a set of grid points (Omega_1).
popGeneratorHS.m: generating true data over a set of grid points (Omega_2).


3. Trivariate Penalized Spline over Triangulation (TPST) estimation:
TPST_est.m: Model fitting via TPST. (GCV is used to choose the best tuning parameter.)


4. Functions for constructing TPST basis, penalty function, and smoothness constraints:
basis.p: generating spline basis.
energyM3D.p: generating energy/penalty function.
smoothness.p: generating H matrix for smoothness conditions.


5. Other functions that helped to construct TPST estimator include:
bary.p, build.p, choose.p, dtri.p, HQbary.p, HQblkdiag.p, indices.p, indices3d.p, insideVT.p, loop3.p, qrH.p, and volume.p.

