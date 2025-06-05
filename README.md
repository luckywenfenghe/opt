MATLAB code for density-based topology optimisation of Navier-Stokes fluid flow with thermal. The code is based on journal publication: "A detailed introduction to density-based topology optimisation of fluid flow problems with implementation in MATLAB", Joe Alexandersen (2023), SMO 66:12, doi:10.1007/s00158-022-03420-9

The authors preprint is available here as "Alexandersen2022_preprint.pdf" and a detailed description of the files of the code base is given in "supplementary_codeDescription.pdf".you can serch it for deep understanding.

The code follow by the https://github.com/sdu-multiphysics/topflow

adaptive mvlim + Heaviside β-projection used for accelerate

Q1: R_T_strong = rho*Cp*(ux(1)*dTdx(1) + ux(2)*dTdx(2)) - Q;
suggest using cuda to accelerate the opt process,Avoid repeatedly creating sparse matrices,
Reuse decomposition / use iterative solution instead :
[Lfac,Ufac,Pfac,Qfac] = lu(J(freedofs,freedofs),'vector');
dS = -Qfac*(Ufac\(Lfac\(Pfac*R(freedofs))));
