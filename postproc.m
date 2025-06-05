%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    NAVIER-STOKES TOPOLOGY OPTIMISATION CODE, MAY 2022    %
% COPYRIGHT (c) 2022, J ALEXANDERSEN. BSD 3-CLAUSE LICENSE %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% POST-PROCESSING
% SETTING UP NODAL COORDINATES FOR STREAMLINES
[X,Y] = meshgrid(0.5:nodx,0.5:nody); sx = 0.5*ones(1,21); sy = linspace(0,nody,21);
U = S(1:(nely+1)*(nelx+1)*2); umag=reshape(sqrt(U(1:2:end).^2+U(2:2:end).^2),nely+1,nelx+1);

% Create filename prefix based on problem type
filename_prefix = sprintf('Problem%d_', probtype);

% DESIGN FIELD
figure(1); imagesc(xPhys); colorbar; caxis([0 1]); axis equal; axis off;
h = streamline(X,Y,reshape(U(1:2:end),nody,nodx),-reshape(U(2:2:end),nody,nodx),sx,sy); set(h,'Color','black');
title('Design Field with Streamlines');
% Save design field
saveas(gcf, [filename_prefix 'design_field.png']);
fprintf('      Saved: %s\n', [filename_prefix 'design_field.png']);

% BRINKMAN PENALTY FACTOR
figure(2); imagesc(reshape(log10(alpha),nely,nelx)); colorbar; caxis([0 log10(alphamax)]); axis equal; axis off; %colormap turbo;
title('Brinkman Penalty Factor (log10)');
% Save Brinkman penalty
saveas(gcf, [filename_prefix 'brinkman_penalty.png']);
fprintf('      Saved: %s\n', [filename_prefix 'brinkman_penalty.png']);

% VELOCITY MAGNITUDE FIELD
figure(3); imagesc(umag); colorbar; axis equal; axis on; hold on; %colormap turbo;
h = streamline(X,Y,reshape(U(1:2:end),nody,nodx),-reshape(U(2:2:end),nody,nodx),sx,sy); set(h,'Color','black');
title('Velocity Magnitude Field with Streamlines');
% Save velocity field
saveas(gcf, [filename_prefix 'velocity_magnitude.png']);
fprintf('      Saved: %s\n', [filename_prefix 'velocity_magnitude.png']);

% PRESSURE FIELD
P = S(2*nodtot+1:3*nodtot);
figure(4); imagesc(reshape(P,nody,nodx)); colorbar; axis equal; axis off; %colormap turbo;
h = streamline(X,Y,reshape(U(1:2:end),nody,nodx),-reshape(U(2:2:end),nody,nodx),sx,sy); set(h,'Color','black');
title('Pressure Field with Streamlines');
% Save pressure field
saveas(gcf, [filename_prefix 'pressure_field.png']);
fprintf('      Saved: %s\n', [filename_prefix 'pressure_field.png']);

% TEMPERATURE FIELD (for problem 3)
if (probtype == 3)
    T = S(3*nodtot+1:4*nodtot);
    Tfield = reshape(T,nody,nodx);
    figure(7); imagesc(Tfield); colorbar; axis equal; axis off; 
    title('Temperature Field (°C)'); colormap('hot');
    h = streamline(X,Y,reshape(U(1:2:end),nody,nodx),-reshape(U(2:2:end),nody,nodx),sx,sy); set(h,'Color','blue');
    % Save temperature field
    saveas(gcf, [filename_prefix 'temperature_field.png']);
    fprintf('      Saved: %s\n', [filename_prefix 'temperature_field.png']);
    
    % Temperature along a line
    if (probtype == 2 || probtype == 3)
        Tline=flipud(diag(fliplr(Tfield))); 
    else
        Tline=Tfield(:,floor((end-1)/2));
    end
    
    % Combined temperature and velocity plots
    figure(8);
    subplot(2,2,1); imagesc(Tfield); colorbar; axis equal; axis off; title('Temperature (°C)'); colormap('hot');
    subplot(2,2,2); imagesc(umag); colorbar; axis equal; axis off; title('Velocity Magnitude'); colormap('jet');
    subplot(2,2,3); plot(Tline,'-r','LineWidth',2); grid on; title('Temperature along diagonal'); xlabel('Position'); ylabel('T (°C)');
    subplot(2,2,4); imagesc(xPhys); colorbar; caxis([0 1]); axis equal; axis off; title('Design Field'); colormap('gray');
    % Save combined temperature/velocity plot
    saveas(gcf, [filename_prefix 'temperature_velocity_combined.png']);
    fprintf('      Saved: %s\n', [filename_prefix 'temperature_velocity_combined.png']);
end

% VELOCITY ALONG A LINE
if (probtype == 2 || probtype == 3)
    uline=flipud(diag(fliplr(umag))); xline=flipud(diag(fliplr(xPhys)));
else
    uline=umag(:,floor((end-1)/2)); xline=xPhys(:,floor(end/2));
end
figure(5);
subplot(3,1,1); plot(uline,'-x'); grid on; title('Velocity magnitude');
subplot(3,1,2); plot(log10(uline),'-x'); grid on; title('Log10(Velocity magnitude)');
subplot(3,1,3); plot(xline,'-x'); grid on; title('Design field'); drawnow
% Save velocity line plots
saveas(gcf, [filename_prefix 'velocity_line_analysis.png']);
fprintf('      Saved: %s\n', [filename_prefix 'velocity_line_analysis.png']);

% TEMPERATURE ITERATION HISTORY (for problem 3)
if (probtype == 3 && loop > 0)
    % Store temperature statistics for iteration history
    if ~exist('Tmax_history','var')
        Tmax_history = [];
        Tmin_history = [];
        Tavg_history = [];
    end
    Tmax_history(end+1) = max(T(:));
    Tmin_history(end+1) = min(T(:));
    Tavg_history(end+1) = mean(T(:));
    
    figure(9);
    subplot(2,1,1); 
    plot(1:length(Tmax_history), Tmax_history, '-r', 'LineWidth', 2); hold on;
    plot(1:length(Tmin_history), Tmin_history, '-b', 'LineWidth', 2);
    plot(1:length(Tavg_history), Tavg_history, '-k', 'LineWidth', 2); hold off;
    legend('T_{max}', 'T_{min}', 'T_{avg}', 'Location', 'best');
    title('Temperature Evolution'); xlabel('Iteration'); ylabel('Temperature (°C)'); grid on;
    
    subplot(2,1,2);
    plot(1:length(Tavg_history), Tavg_history - 25, '-g', 'LineWidth', 2);
    title('Temperature Rise from Inlet'); xlabel('Iteration'); ylabel('ΔT (°C)'); grid on;
    % Save temperature evolution plot
    saveas(gcf, [filename_prefix 'temperature_evolution.png']);
    fprintf('      Saved: %s\n', [filename_prefix 'temperature_evolution.png']);
end

%% SUMMARY OF SAVED FILES
fprintf('=========================================================\n');
fprintf('      POST-PROCESSING COMPLETE - FILES SAVED:\n');
fprintf('=========================================================\n');
if (probtype == 3)
    fprintf('      Total files saved: 7 images\n');
    fprintf('      • Design field with streamlines\n');
    fprintf('      • Brinkman penalty factor\n');
    fprintf('      • Velocity magnitude field\n');
    fprintf('      • Pressure field\n');
    fprintf('      • Velocity line analysis\n');
    fprintf('      • Temperature field\n');
    fprintf('      • Temperature/velocity combined view\n');
    if (loop > 0)
        fprintf('      • Temperature evolution history\n');
    end
else
    fprintf('      Total files saved: 5 images\n');
    fprintf('      • Design field with streamlines\n');
    fprintf('      • Brinkman penalty factor\n');
    fprintf('      • Velocity magnitude field\n');
    fprintf('      • Pressure field\n');
    fprintf('      • Velocity line analysis\n');
end
fprintf('      All files saved with prefix: %s\n', filename_prefix);
fprintf('=========================================================\n');
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