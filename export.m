%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    NAVIER-STOKES TOPOLOGY OPTIMISATION CODE, MAY 2022    %
% COPYRIGHT (c) 2022, J ALEXANDERSEN. BSD 3-CLAUSE LICENSE %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DXF
xMod = padarray(padarray(xPhys,[1 1],'replicate'),[1 1],1);
figure('visible','off');
C = contour(flipud(xMod),[0.5 0.5]);
C = unique(C','rows','stable')';
C = C-2.5; C(1,:) = C(1,:)*dx; C(2,:) = C(2,:)*dy;
h = sqrt(dx^2+dy^2);
ind = 1:size(C,2);
ind1 = intersect(find(C(1,:)<Lx+2*dx),find(C(2,:)<Ly+2*dy));
ind2 = intersect(find(C(1,:)>0-dx),find(C(2,:)>0-dy));
ind = intersect(ind1,ind2);
X = []; Y = []; c = 0;
indused = [];
for i = 1:(length(ind))
    x0 = C(1,ind(i)); y0 = C(2,ind(i));
    redind = ind([1:max(1,i-5) i+1:min(length(ind),i+2)]);
    redind = setdiff(redind,indused(max(1,c-10):end));
    R = sqrt( (C(1,redind)-x0).^2 + (C(2,redind)-y0).^2 );
    [minR,idx] = min(R);
    if (minR < 2*h && minR > 0)
        c = c + 1;
        indused(c) = i;
        X(c,:) = [x0 C(1,redind(idx))];
        Y(c,:) = [y0 C(2,redind(idx))];
    end
end
fid=fopen([filename '.dxf'],'w');
fprintf(fid,'0\nSECTION\n2\nENTITIES\n0\n');
for i=1:size(X,1)
      fprintf(fid,'LINE\n8\n0\n');
      fprintf(fid,'10\n%.4f\n20\n%.4f\n30\n%.4f\n',X(i,1),Y(i,1),0);
      fprintf(fid,'11\n%.4f\n21\n%.4f\n31\n%.4f\n',X(i,2),Y(i,2),0); 
      fprintf(fid,'0\n');
end
fprintf(fid,'ENDSEC\n0\nEOF\n');
fclose(fid);

%% EXPORT SOLUTION DATA (CSV format)
% Export velocity field
U = S(1:(nely+1)*(nelx+1)*2); 
Ux = reshape(U(1:2:end),nely+1,nelx+1);
Uy = reshape(U(2:2:end),nely+1,nelx+1);
umag = sqrt(Ux.^2 + Uy.^2);

% Export pressure field
P = S(2*nodtot+1:3*nodtot);
Pfield = reshape(P,nody,nodx);

% Create coordinate matrices
[Xgrid,Ygrid] = meshgrid(linspace(0,Lx,nelx+1),linspace(0,Ly,nely+1));

% Export velocity and pressure data
fid = fopen([filename '_flow.csv'],'w');
fprintf(fid,'X,Y,Ux,Uy,Umag,Pressure,Design\n');
for i = 1:size(Xgrid,1)
    for j = 1:size(Xgrid,2)
        if i <= size(Pfield,1) && j <= size(Pfield,2)
            design_val = 0; % Default for boundary nodes
            if i <= nely && j <= nelx
                design_val = xPhys(i,j);
            end
            fprintf(fid,'%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n',...
                    Xgrid(i,j),Ygrid(i,j),Ux(i,j),Uy(i,j),umag(i,j),Pfield(i,j),design_val);
        end
    end
end
fclose(fid);

% Export temperature field data (for problem 3)
if (probtype == 3)
    T = S(3*nodtot+1:4*nodtot);
    Tfield = reshape(T,nody,nodx);
    
    fid = fopen([filename '_temperature.csv'],'w');
    fprintf(fid,'X,Y,Temperature,Design\n');
    for i = 1:size(Xgrid,1)
        for j = 1:size(Xgrid,2)
            if i <= size(Tfield,1) && j <= size(Tfield,2)
                design_val = 0; % Default for boundary nodes
                if i <= nely && j <= nelx
                    design_val = xPhys(i,j);
                end
                fprintf(fid,'%.6f,%.6f,%.6f,%.6f\n',...
                        Xgrid(i,j),Ygrid(i,j),Tfield(i,j),design_val);
            end
        end
    end
    fclose(fid);
    
    % Export temperature statistics
    fid = fopen([filename '_temp_stats.txt'],'w');
    fprintf(fid,'Temperature Field Statistics\n');
    fprintf(fid,'============================\n');
    fprintf(fid,'Minimum Temperature: %.2f °C\n', min(T(:)));
    fprintf(fid,'Maximum Temperature: %.2f °C\n', max(T(:)));
    fprintf(fid,'Average Temperature: %.2f °C\n', mean(T(:)));
    fprintf(fid,'Temperature Rise: %.2f °C\n', mean(T(:)) - 25);
    fprintf(fid,'Heat Source: %.1f W/m³\n', Qsource);
    fprintf(fid,'Thermal Conductivity: %.2f W/m·K\n', k);
    fprintf(fid,'Specific Heat: %.0f J/kg·K\n', Cp);
    fclose(fid);
end

%% EXPORT MATLAB DATA
save([filename '_results.mat'], 'xPhys', 'S', 'obj', 'loop', 'Renum');
if (probtype == 3)
    save([filename '_results.mat'], 'xPhys', 'S', 'obj', 'loop', 'Renum', 'Tfield', 'k', 'Cp', 'Qsource', '-append');
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code was initial written by: Joe Alexandersen                      %
%                           Department of Mechanical and                  %
%                                         Electrical Engineering          %
%                           University of Southern Denmark                %
%                           DK-5230 Odense M, Denmark.                    % 
% Has been refined by authors: https://github.com/luckywenfenghe          %
% Adaptive Move-Limit via Trust-Region MMA                                %
% Automatic β-Projection Continuation                                     %
% Re-use of Sparse LU Factorization in Newton Loops                       %
%                                                                         %
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