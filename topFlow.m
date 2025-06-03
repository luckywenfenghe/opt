%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    NAVIER-STOKES TOPOLOGY OPTIMISATION CODE, MAY 2022    %
% COPYRIGHT (c) 2022, J ALEXANDERSEN. BSD 3-CLAUSE LICENSE %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc;
%% DEFINITION OF INPUT PARAMETERS
% PROBLEM TO SOLVE (1 = DOUBLE PIPE; 2 = PIPE BEND; 3 = PIPE BEND WITH HEAT TRANSFER)
probtype = 3;
% DOMAIN SIZE
Lx = 1.0; Ly = 1.0;
% DOMAIN DISCRETISATION
nely = 40; nelx = nely*Lx/Ly;
% ALLOWABLE FLUID VOLUME FRACTION
volfrac = 1/4; xinit = volfrac;
% PHYSICAL PARAMETERS
Uin = 1e1; rho = 1e1; mu = 1e1;
% THERMAL PARAMETERS (for problem 3)
if (probtype == 3)
    kappa = 0.4; % Thermal conductivity (W/m·K) - renamed from k to avoid conflicts
    Cp = 4180; % Specific heat capacity (J/kg·K)
    dt_thermal = 0.01; % Time step for transient thermal analysis
else
    kappa = 1; Cp = 1; dt_thermal = 1; % Default values for non-thermal problems
end
% BRINKMAN PENALISATION
alphamax = 2.5*mu/(0.01^2); alphamin = 2.5*mu/(100^2);
% CONTINUATION STRATEGY
ainit = 2.5*mu/(0.1^2);
qinit = (-xinit*(alphamax-alphamin) - ainit + alphamax)/(xinit*(ainit-alphamin));
qavec = qinit./[1 2 10 20]; qanum = length(qavec); conit = 50;
% OPTIMISATION PARAMETERS
maxiter = qanum*conit; mvlim = 0.2; plotdes = 1;
chlim = 1e-3; chnum = 5;
% NEWTON SOLVER PARAMETERS
nltol = 1e-6; nlmax = 25; plotres = 0;
% EXPORT FILE
filename='output'; exportdxf = 0;
%% PREPARE FINITE ELEMENT ANALYSIS
dx = Lx/nelx; dy = Ly/nely;
nodx = nelx+1; nody = nely+1; nodtot = nodx*nody;
neltot = nelx*nely; 
if (probtype == 3)
    doftot = 4*nodtot; % u(2) + p(1) + T(1) = 4 DOFs per node
else
    doftot = 3*nodtot; % u(2) + p(1) = 3 DOFs per node
end
% NODAL CONNECTIVITY
nodenrs = reshape(1:nodtot,nody,nodx);
edofVecU = reshape(2*nodenrs(1:end-1,1:end-1)+1,neltot,1);
edofMatU = repmat(edofVecU,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],neltot,1);
edofVecP = reshape(nodenrs(1:end-1,1:end-1),neltot,1);
edofMatP = repmat(edofVecP,1,4)+repmat([1 nely+[2 1] 0],neltot,1);
if (probtype == 3)
    edofVecT = reshape(nodenrs(1:end-1,1:end-1),neltot,1);
    edofMatT = repmat(edofVecT,1,4)+repmat([1 nely+[2 1] 0],neltot,1);
    edofMat = [edofMatU 2*nodtot+edofMatP 3*nodtot+edofMatT];
    iJ = reshape(kron(edofMat,ones(16,1))',256*neltot,1);
    jJ = reshape(kron(edofMat,ones(1,16))',256*neltot,1);
    iR = reshape(edofMat',16*neltot,1); jR = ones(16*neltot,1); 
    jE = repmat(1:neltot,16,1);
else
    edofMat = [edofMatU 2*nodtot+edofMatP];
    iJ = reshape(kron(edofMat,ones(12,1))',144*neltot,1);
    jJ = reshape(kron(edofMat,ones(1,12))',144*neltot,1);
    iR = reshape(edofMat',12*neltot,1); jR = ones(12*neltot,1); 
    jE = repmat(1:neltot,12,1);
end
%% DEFINE BOUNDARY CONDITIONS
% DEFINE THE PROBLEMS IN SEPARATE MATLAB FILE
run('problems.m');
% NULLSPACE MATRICES FOR IMPOSING BOUNDARY CONDITIONS
EN=speye(doftot); ND=EN; ND(fixedDofs,fixedDofs)=0.0; EN=EN-ND;
% VECTORS FOR FREE DOFS
alldofs = 1:doftot; freedofs = setdiff(alldofs,fixedDofs);
%% INITIALISATION
% SOLUTION VECTOR
S = zeros(doftot,1); dS = S; L = S; 
S(fixedDofs) = DIR(fixedDofs);
% PREVIOUS TEMPERATURE FIELD (for transient analysis)
if (probtype == 3)
    Sold = S; % Initialize previous solution
end
% DESIGN FIELD
xPhys = xinit*ones(nely,nelx); 
% COUNTERS
loop = 0; loopcont = 0; nlittot = 0; chcnt = 0;
% CHANGE
change = Inf; objOld = Inf;
% CONTINUATION
qastep = 1; qa = qavec(1);
% VECTORISED CONSTANTS
dxv = dx*ones(1,neltot); dyv = dy*ones(1,neltot);
muv = mu*ones(1,neltot); rhov = rho*ones(1,neltot);
dtv = dt_thermal*ones(1,neltot); % Time step vector
if (probtype == 3)
    kv = kappa*ones(1,neltot); Cpv = Cp*ones(1,neltot);
    Qv = Qsource*ones(1,neltot);
end
%% OUTPUT PROBLEM INFORMATION
fprintf('=========================================================\n');
fprintf('      Problem number: %2i - Reynolds number: %3.2e\n',probtype,Renum);
if (probtype == 3)
    fprintf('      Heat transfer problem with uniform heat source\n');
end
fprintf('=========================================================\n');
fprintf('      Design it.:   0\n');
%% START ITERATION
destime = tic; ittime = tic;
while (loop <= maxiter)
    if (plotdes); figure(1); imagesc(xPhys); colorbar; caxis([0 1]); axis equal; axis off; drawnow; end
    %% GREYSCALE INDICATOR
    Md = 100*full(4*sum(xPhys(:).*(1-xPhys(:)))/neltot); 
    %% MATERIAL INTERPOLATION
    alpha = alphamin + (alphamax-alphamin)*(1-xPhys(:))./(1+qa*xPhys(:));
    dalpha = (qa*(alphamax - alphamin)*(xPhys(:) - 1))./(xPhys(:)*qa + 1).^2 - (alphamax - alphamin)./(xPhys(:)*qa + 1);
    %% NON-LINEAR NEWTON SOLVER
    normR = 1; nlit = 0; fail = -1; nltime = tic;
    while (fail ~= 1)
        nlit = nlit+1; nlittot = nlittot+1;
        % BUILD RESIDUAL AND JACOBIAN
        if (probtype == 3)
            % Extract state variables for thermal problem (16 DOFs per element)
            uVars = S(edofMat(:,1:8)); % Velocity components u1-u8
            pVars = S(edofMat(:,9:12)); % Pressure components p1-p4  
            TVars = S(edofMat(:,13:16)); % Temperature components T1-T4
            
            sR = zeros(16, neltot);
            for i = 1:neltot
                sR(:,i) = RES(dxv(i),dyv(i),muv(i),rhov(i),kv(i),Cpv(i),Qv(i),alpha(i),...
                    uVars(i,1),uVars(i,2),uVars(i,3),uVars(i,4),uVars(i,5),uVars(i,6),uVars(i,7),uVars(i,8),...
                    pVars(i,1),pVars(i,2),pVars(i,3),pVars(i,4),...
                    TVars(i,1),TVars(i,2),TVars(i,3),TVars(i,4));
            end
        else
            % Extract state variables for non-thermal problem (12 DOFs per element)
            uVars = S(edofMat(:,1:8)); % Velocity components u1-u8
            pVars = S(edofMat(:,9:12)); % Pressure components p1-p4
            
            sR = zeros(16, neltot); % Still need 16 DOFs for consistent function interface
            for i = 1:neltot
                sR(:,i) = RES(dxv(i),dyv(i),muv(i),rhov(i),kv(i),Cpv(i),Qv(i),alpha(i),...
                    uVars(i,1),uVars(i,2),uVars(i,3),uVars(i,4),uVars(i,5),uVars(i,6),uVars(i,7),uVars(i,8),...
                    pVars(i,1),pVars(i,2),pVars(i,3),pVars(i,4),...
                    0,0,0,0); % Default temperature values for non-thermal problems
            end
            sR = sR(1:12,:); % Extract only momentum and continuity equations
        end
        R = sparse(iR,jR,sR(:)); R(fixedDofs) = 0;
        if (nlit == 1); r0 = norm(R); end
        r1 = norm(R); normR = r1/r0;
        if (plotres); figure(6); semilogy(nlittot,normR,'x'); axis square; grid on; hold on; 
        end
        if (normR < nltol); break;
        end
        
        if (probtype == 3)
            sJ = zeros(256, neltot);
            for i = 1:neltot
                sJ(:,i) = JAC(dxv(i),dyv(i),muv(i),rhov(i),kv(i),Cpv(i),Qv(i),alpha(i),...
                    uVars(i,1),uVars(i,2),uVars(i,3),uVars(i,4),uVars(i,5),uVars(i,6),uVars(i,7),uVars(i,8),...
                    pVars(i,1),pVars(i,2),pVars(i,3),pVars(i,4),...
                    TVars(i,1),TVars(i,2),TVars(i,3),TVars(i,4));
                end

        else
            sJ = zeros(256, neltot); % Still need 256 for consistent function interface
            for i = 1:neltot
                sJ(:,i) = JAC(dxv(i),dyv(i),muv(i),rhov(i),kv(i),Cpv(i),Qv(i),alpha(i),...
                    uVars(i,1),uVars(i,2),uVars(i,3),uVars(i,4),uVars(i,5),uVars(i,6),uVars(i,7),uVars(i,8),...
                    pVars(i,1),pVars(i,2),pVars(i,3),pVars(i,4),...
                    0,0,0,0); % Default temperature values for non-thermal problems
                end

                sJ = sJ(1:144,:); % Extract only momentum and continuity entries
        end
        J = sparse(iJ,jJ,sJ(:)); J = (ND'*J*ND+EN);  
        % CALCULATE NEWTON STEP
        dS = -J\R;
        % L2-NORM LINE SEARCH
        Sp = S + 0.5*dS;
        if (probtype == 3)
            uVars = Sp(edofMat(:,1:8));
            pVars = Sp(edofMat(:,9:12));
            TVars = Sp(edofMat(:,13:16));

            sR = zeros(16, neltot);
            for i = 1:neltot
                sR(:,i) = RES(dxv(i),dyv(i),muv(i),rhov(i),kv(i),Cpv(i),Qv(i),alpha(i),...
                    uVars(i,1),uVars(i,2),uVars(i,3),uVars(i,4),uVars(i,5),uVars(i,6),uVars(i,7),uVars(i,8),...
                    pVars(i,1),pVars(i,2),pVars(i,3),pVars(i,4),...
                    TVars(i,1),TVars(i,2),TVars(i,3),TVars(i,4));
                end
            
        else
            uVars = Sp(edofMat(:,1:8));
            pVars = Sp(edofMat(:,9:12));

            sR = zeros(12, neltot);
            for i = 1:neltot
                sR(:,i) = RES(dxv(i),dyv(i),muv(i),rhov(i),alpha(i),...
                    uVars(i,1),uVars(i,2),uVars(i,3),uVars(i,4),uVars(i,5),uVars(i,6),uVars(i,7),uVars(i,8),...
                    pVars(i,1),pVars(i,2),pVars(i,3),pVars(i,4));
                end
            end
        
        R = sparse(iR,jR,sR(:)); R(fixedDofs) = 0; r2 = norm(R);
        Sp = S + 1.0*dS;
        if (probtype == 3)
            uVars = Sp(edofMat(:,1:8));
            pVars = Sp(edofMat(:,9:12));
            TVars = Sp(edofMat(:,13:16));

            sR = zeros(16, neltot);
            for i = 1:neltot
                sR(:,i) = RES(dxv(i),dyv(i),muv(i),rhov(i),kv(i),Cpv(i),Qv(i),alpha(i),...
                    uVars(i,1),uVars(i,2),uVars(i,3),uVars(i,4),uVars(i,5),uVars(i,6),uVars(i,7),uVars(i,8),...
                    pVars(i,1),pVars(i,2),pVars(i,3),pVars(i,4),...
                    TVars(i,1),TVars(i,2),TVars(i,3),TVars(i,4));
                end
            
        else
            uVars = Sp(edofMat(:,1:8));
            pVars = Sp(edofMat(:,9:12));

            sR = zeros(12, neltot);
            for i = 1:neltot
                sR(:,i) = RES(dxv(i),dyv(i),muv(i),rhov(i),alpha(i),...
                    uVars(i,1),uVars(i,2),uVars(i,3),uVars(i,4),uVars(i,5),uVars(i,6),uVars(i,7),uVars(i,8),...
                    pVars(i,1),pVars(i,2),pVars(i,3),pVars(i,4));
                end
            end
        
        R = sparse(iR,jR,sR(:)); R(fixedDofs) = 0; r3 = norm(R);
        % SOLUTION UPDATE WITH "OPTIMAL" DAMPING
        lambda = max(0.01,min(1.0,(3*r1 + r3 - 4*r2)/(4*r1 + 4*r3 - 8*r2)));
        S = S + lambda*dS;
        % IF FAIL, RETRY FROM ZERO SOLUTION
        if (nlit == nlmax && fail < 0); nlit = 0; S(freedofs) = 0.0; normR=1; fail = fail+1; end
        if (nlit == nlmax && fail < 1); fail = fail+1; end
    end
    nltime=toc(nltime);
    fprintf('      Newton it.: %2i - Res. norm: %3.2e - Sol. time: %6.3f sec\n',nlit,normR,nltime);
    if (fail == 1); error('ERROR: Solver did not converge after retry from zero!\n      Stopping optimisation.\n'); end
    %% OBJECTIVE EVALUATION
    if (probtype == 3)
        phiVals = zeros(1, neltot);
        for i = 1:neltot
            phiVals(i) = PHI(dxv(i),dyv(i),muv(i),kv(i),alpha(i),...
                uVars(i,1),uVars(i,2),uVars(i,3),uVars(i,4),uVars(i,5),uVars(i,6),uVars(i,7),uVars(i,8),...
                pVars(i,1),pVars(i,2),pVars(i,3),pVars(i,4),...
                TVars(i,1),TVars(i,2),TVars(i,3),TVars(i,4));
        end
        obj = sum(phiVals);
    else
        phiVals = zeros(1, neltot);
        for i = 1:neltot
            phiVals(i) = PHI(dxv(i),dyv(i),muv(i),alpha(i),...
                uVars(i,1),uVars(i,2),uVars(i,3),uVars(i,4),uVars(i,5),uVars(i,6),uVars(i,7),uVars(i,8),...
                pVars(i,1),pVars(i,2),pVars(i,3),pVars(i,4));
        end
        obj = sum(phiVals);
    end
    change = abs(objOld-obj)/objOld; objOld = obj;
    %% VOLUME CONSTRAINT
    V = mean(xPhys(:));
    %% PRINT RESULTS
    ittime = toc(ittime);
    fprintf('      Obj.: %3.2e - Constr.: %3.2e - Md: %3.2f\n',obj,V,Md);
    fprintf('      Change: %4.3e - It. time: %6.3f sec\n',change,ittime);
    fprintf('      Contin. step: %2i - qa: %4.3e\n',qastep,qa);
    ittime = tic;
    %% EVALUATE CURRENT ITERATE - CONTINUE UNLESS CONSIDERED CONVERGED
    if (change < chlim); chcnt = chcnt + 1; else; chcnt = 0; end
    if (qastep == qanum && ( (chcnt == chnum) || (loopcont == conit) ) ); end
    %% PRINT HEADER FOR ITERATION
    loop = loop + 1; loopcont = loopcont + 1;
    fprintf('---------------------------------------------------------\n');
    fprintf('      Design it.:%4i\n',loop);
    %% ADJOINT SOLVER
    if (probtype == 3)
        dphiVals = zeros(16, neltot);
        for i = 1:neltot
            dphiVals(:,i) = dPHIds(dxv(i),dyv(i),muv(i),kv(i),alpha(i),...
                uVars(i,1),uVars(i,2),uVars(i,3),uVars(i,4),uVars(i,5),uVars(i,6),uVars(i,7),uVars(i,8),...
                pVars(i,1),pVars(i,2),pVars(i,3),pVars(i,4),...
                TVars(i,1),TVars(i,2),TVars(i,3),TVars(i,4));
        end
        sR = [dphiVals; zeros(0,neltot)]; % No padding needed as dPHIds now returns 16 components
    else
        dphiVals = zeros(12, neltot);
        for i = 1:neltot
            dphiVals(:,i) = dPHIds(dxv(i),dyv(i),muv(i),alpha(i),...
                uVars(i,1),uVars(i,2),uVars(i,3),uVars(i,4),uVars(i,5),uVars(i,6),uVars(i,7),uVars(i,8),...
                pVars(i,1),pVars(i,2),pVars(i,3),pVars(i,4));
        end
        sR = [dphiVals; zeros(0,neltot)]; % No padding needed
    end
    RHS = sparse(iR,jR,sR(:)); RHS(fixedDofs) = 0;
    if (probtype == 3)
        sJ = zeros(256, neltot);
        for i = 1:neltot
            sJ(:,i) = JAC(dxv(i),dyv(i),muv(i),rhov(i),kv(i),Cpv(i),Qv(i),alpha(i),...
                uVars(i,1),uVars(i,2),uVars(i,3),uVars(i,4),uVars(i,5),uVars(i,6),uVars(i,7),uVars(i,8),...
                pVars(i,1),pVars(i,2),pVars(i,3),pVars(i,4),...
                TVars(i,1),TVars(i,2),TVars(i,3),TVars(i,4));
        end
    else
        sJ = zeros(256, neltot); % Still need 256 for consistent function interface
        for i = 1:neltot
            sJ(:,i) = JAC(dxv(i),dyv(i),muv(i),rhov(i),kv(i),Cpv(i),Qv(i),alpha(i),...
                uVars(i,1),uVars(i,2),uVars(i,3),uVars(i,4),uVars(i,5),uVars(i,6),uVars(i,7),uVars(i,8),...
                pVars(i,1),pVars(i,2),pVars(i,3),pVars(i,4),...
                0,0,0,0); % Default temperature values for non-thermal problems
        end
        sJ = sJ(1:144,:); % Extract only momentum and continuity entries
    end
    J = sparse(iJ,jJ,sJ(:)); J = (ND'*J*ND+EN);
    L = J'\RHS;
    %% COMPUTE SENSITIVITIES
    % OBJECTIVE
    if (probtype == 3)
        drdgVals = zeros(16, neltot);
        dphidgVals = zeros(1, neltot);
        for i = 1:neltot
            drdgVals(:,i) = dRESdg(dxv(i),dyv(i),muv(i),rhov(i),kv(i),Cpv(i),Qv(i),alpha(i),dalpha(i),...
                uVars(i,1),uVars(i,2),uVars(i,3),uVars(i,4),uVars(i,5),uVars(i,6),uVars(i,7),uVars(i,8),...
                pVars(i,1),pVars(i,2),pVars(i,3),pVars(i,4),...
                TVars(i,1),TVars(i,2),TVars(i,3),TVars(i,4));
            dphidgVals(i) = dPHIdg(dxv(i),dyv(i),muv(i),kv(i),alpha(i),dalpha(i),...
                uVars(i,1),uVars(i,2),uVars(i,3),uVars(i,4),uVars(i,5),uVars(i,6),uVars(i,7),uVars(i,8),...
                pVars(i,1),pVars(i,2),pVars(i,3),pVars(i,4),...
                TVars(i,1),TVars(i,2),TVars(i,3),TVars(i,4));
        end
        dRdg = sparse(iR(:),jE(:),drdgVals(:));
        dphidg = dphidgVals;
        % Add temperature equation contribution to sensitivity
        try
            sR_temp = zeros(4, neltot);
            for i = 1:neltot
                sR_temp(:,i) = dRESTemperature(dxv(i),dyv(i),muv(i),rhov(i),kv(i),Cpv(i),Qv(i),alpha(i),...
                    uVars(i,1),uVars(i,2),uVars(i,3),uVars(i,4),uVars(i,5),uVars(i,6),uVars(i,7),uVars(i,8),...
                    pVars(i,1),pVars(i,2),pVars(i,3),pVars(i,4),...
                    TVars(i,1),TVars(i,2),TVars(i,3),TVars(i,4));
            end
            % Extract temperature adjoint variables (last 4 entries per element)
            L_temp_elem = L(edofMat(:,13:16)); % Extract 4 temperature adjoint values per element (neltot×4)
            L_temp = L_temp_elem.'; % Transpose to 4×neltot to match sR_temp dimensions
            temp_contrib = sum(L_temp .* sR_temp, 1); % Calculate Σ λ_T · ∂R_T/∂ρ per element (1 × neltot)
            % Add temperature contribution to design sensitivity
            sens = reshape(dphidg - L'*dRdg,nely,nelx) + reshape(temp_contrib, nely, nelx);
        catch ME
            % If dRESTemperature doesn't exist yet, use standard sensitivity
            sens = reshape(dphidg - L'*dRdg,nely,nelx);
            fprintf('      Warning: dRESTemperature.m call failed (%s), using standard sensitivity\n', ME.message);
        end
    else
        drdgVals = zeros(12, neltot);
        dphidgVals = zeros(1, neltot);
        for i = 1:neltot
            drdgVals(:,i) = dRESdg(dxv(i),dyv(i),muv(i),rhov(i),alpha(i),dalpha(i),...
                uVars(i,1),uVars(i,2),uVars(i,3),uVars(i,4),uVars(i,5),uVars(i,6),uVars(i,7),uVars(i,8),...
                pVars(i,1),pVars(i,2),pVars(i,3),pVars(i,4));
            dphidgVals(i) = dPHIdg(dxv(i),dyv(i),muv(i),alpha(i),dalpha(i),...
                uVars(i,1),uVars(i,2),uVars(i,3),uVars(i,4),uVars(i,5),uVars(i,6),uVars(i,7),uVars(i,8),...
                pVars(i,1),pVars(i,2),pVars(i,3),pVars(i,4));
        end
        dRdg = sparse(iR(:),jE(:),drdgVals(:));
        dphidg = dphidgVals;
        sens = reshape(dphidg - L'*dRdg,nely,nelx);
    end
    % VOLUME CONSTRAINT
    dV = ones(nely,nelx)/neltot;
    %% OPTIMALITY CRITERIA UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES
    xnew = xPhys; xlow = xPhys(:)-mvlim; xupp = xPhys(:)+mvlim;
    ocfac = xPhys(:).*max(1e-10,(-sens(:)./dV(:))).^(1/3);
    l1 = 0; l2 = ( 1/(neltot*volfrac)*sum( ocfac ) )^3;
    while (l2-l1)/(l1+l2) > 1e-3
        lmid = 0.5*(l2+l1);
        xnew(:) = max(0,max(xlow,min(1,min(xupp,ocfac/(lmid^(1/3))))));
        if mean(xnew(:)) > volfrac; l1 = lmid; else; l2 = lmid; end
    end
    xPhys = xnew;
    %% CONTINUATION UPDATE
    if (qastep < qanum && (loopcont == conit || chcnt == chnum) )
        loopcont = 0; chcnt = 0;
        qastep = qastep + 1; qa = qavec(qastep);
    end
end
%% PRINT FINAL INFORMATION
destime = toc(destime);
fprintf('=========================================================\n');
fprintf('      Number of design iterations: %4i\n',loop);
fprintf('      Final objective: %4.3e\n',obj);
fprintf('      Total time taken: %6.2f min\n',destime/60);
fprintf('=========================================================\n');
%% PLOT RESULTS
run('postproc.m');
if (exportdxf); run('export.m'); end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code was written by: Joe Alexandersen                              %
%                           Department of Mechanical and                  %
%                                         Electrical Engineering          %
%                           University of Southern Denmark                %
%                           DK-5230 Odense M, Denmark.                    %
% Please send your comments and questions to: joal@sdu.dk                 %
%                                                                         %
% The code is intended for educational purposes and theoretical details   %
% are discussed in the paper: "A detailed introduction to density-based   %
% topology optimisation of fluid flow problems including implementation   %
% in MATLAB", J. Alexandersen, SMO 2022, doi:                             %                          
%                                                                         %
% A preprint version of the paper can be downloaded from the author's     %
% website: joealexandersen.com                                            %
% The code is available from GitHub: github.com/sdu-multiphysics/topflow  %
%                                                                         %
% The basic structure of the code is based on the 88-line code for        %
% elastic compliance from: "Efficient topology optimization in MATLAB     %
% using 88 lines of code", E. Andreassen, A. Clausen, M. Schevenels,      %
% B. S. Lazarov and O. Sigmund, SMO 2010, doi:10.1007/s00158-010-0594-7   %
%                                                                         %
% Disclaimer:                                                             %
% The author does not guarantee that the code is free from errors.        %
% Furthermore, the author shall not be liable in any event caused by the  %
% use of the program.                                                     %      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

