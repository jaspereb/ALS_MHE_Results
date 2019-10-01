% This script loads the ALS estimates for Q,R. Then runs the original EKF,
% the EKF-ALS-EKF and the MHE-ALS-EKF to plot all three of these together
% against the ground truth. This script also generates the error tables.

% Innovation frequency plots for all filters are included.

% Note that this code was thrown together from various scripts to generate 
% the paper figures in one go, it is not clean or thoroughly checked. 

% For any discrepancies with the method reported in the paper, you should 
% trust the paper version.

clearvars;
close all;
load('./ALS_Results.mat');

%% Redefine REST or QEST
REST = 1e-3*eye(p); %This comes from the grid search over R,Q
%Qest should come from ALS with Me=21, N=20 MHE data

%% MHE-ALS-EKF
%Define system matrices
functX = @RealXFunction;
functH = @RealGFunction;
R{1} = REST;
Q{1} = QEST;
P{1} = 1.0*eye(size(x,1));

%The noise contribution matrices
G = eye(12);
H = eye(10);

%Initial state
xHatEKF{1} = [0.08 0 0.03 0 0.32 0 0 0 0 0 -pi() 0]';
K = [];
Y = [];
[z,C{1}] = calcJac(functH, xHatEKF{1});

for time = 2:size(y,2)
    %Predict Stage
    [xHatEKF{time},A] = calcJac(functX, xHatEKF{time-1});
    xHatEKFPrior{time} = xHatEKF{time};
    P{time} = A*P{time-1}*A' + Q{time-1};
    
    %Keep R and Q
    R{time} = R{time-1};
    Q{time} = Q{time-1};
    
    [z,C{time}] = calcJac(functH, xHatEKF{time});
    K{time} = P{time}*C{time}'*inv(R{time} + C{time}*P{time}*C{time}');
    xHatEKF{time} = xHatEKF{time} + K{time}*(y(:,time) - z);
    Y_MHE_ALS_EKF{time} = y(:,time) - z;
    P{time} = P{time} - K{time}*C{time}*P{time};
end

MHE_ALS_EKF = xHatEKFPrior;

%MHE-ALS-EKF Error
[meanEKFErrors, maxEKFErrors, stdEKFErrors] = calcEKFError(MHE_ALS_EKF, x);
fprintf('Calculated errors for MHE-ALS-EKF filter: \n');
fprintf('mean: %e  \n', meanEKFErrors);
fprintf('max: %e \n', maxEKFErrors);
fprintf('std: %e \n', stdEKFErrors);

xGT = x;
yGT = y;

clearvars -except xGT yGT MHE_ALS_EKF Y_MHE_ALS_EKF

%% Run the original EKF
measurements = yGT;
rVal = 1e-3;
Qblock = diag([1e-7, 1e-7]);

%Define system matrices
functX = @RealXFunction;
functH = @RealGFunction;

R{1} = rVal*eye(10);
Q{1} = blkdiag(Qblock,Qblock,Qblock,Qblock,Qblock,Qblock);
H = eye(10);

%Initial state
x{1} = [0.08 0 0.03 0 0.32 0 0 0 0 0 -pi() 0]';
P{1} = 1.0*eye(size(x,1));
K = [];
Y = [];
[z,C{1}] = calcJac(functH, x{1});

for time = 2:size(measurements,2)
    %Predict Stage
    [x{time},A] = calcJac(functX, x{time-1});
    xPrior{time} = x{time};
    P{time} = A*P{time-1}*A' + Q{time-1};
    
    %Keep R and Q
    R{time} = R{time-1};
    Q{time} = Q{time-1};
    
    [z,C{time}] = calcJac(functH, x{time});
    K{time} = P{time}*C{time}'*inv(R{time} + C{time}*P{time}*C{time}');
    x{time} = x{time} + K{time}*(measurements(:,time) - z);
    Y_Original_EKF{time} = yGT(:,time) - z;
    P{time} = P{time} - K{time}*C{time}*P{time};
    
    %Precalculate for ALS
    Y{time} = measurements(:,time) - z;
    AK{time} = -A*K{time};
    A_bar{time}=A+AK{time}*C{time};
    AK{time} = AK{time}*H;
end

Original_EKF = xPrior;

%MHE-ALS-EKF Error
[meanEKFErrors, maxEKFErrors, stdEKFErrors] = calcEKFError(xPrior, xGT);
fprintf('Calculated errors for original-EKF filter: \n');
fprintf('mean: %e  \n', meanEKFErrors);
fprintf('max: %e \n', maxEKFErrors);
fprintf('std: %e \n', stdEKFErrors);

clearvars -except xGT yGT MHE_ALS_EKF Y_MHE_ALS_EKF Original_EKF Y_Original_EKF

%% EKF-ALS-EKF
y = yGT;
rVal = 1e-3;
Q{1} = (10^-5)*diag([0.044, 4.0659, 0.044, 4.059, 0.044, 4.059, 1.368, 71.036, 1.368, 71.036, 1.368, 71.036]);
R{1} = rVal*eye(10);

%Define system matrices
functX = @RealXFunction;
functH = @RealGFunction;
P{1} = 1.0*eye(size(xGT,1));

%The noise contribution matrices
G = eye(12);
H = eye(10);

%Initial state
xHatEKF{1} = [0.08 0 0.03 0 0.32 0 0 0 0 0 -pi() 0]';
K = [];
Y = [];
[z,C{1}] = calcJac(functH, xHatEKF{1});

for time = 2:size(y,2)
    %Predict Stage
    [xHatEKF{time},A] = calcJac(functX, xHatEKF{time-1});
    xHatEKFPrior{time} = xHatEKF{time};
    P{time} = A*P{time-1}*A' + Q{time-1};
    
    %Keep R and Q
    R{time} = R{time-1};
    Q{time} = Q{time-1};
    
    [z,C{time}] = calcJac(functH, xHatEKF{time});
    K{time} = P{time}*C{time}'*inv(R{time} + C{time}*P{time}*C{time}');
    xHatEKF{time} = xHatEKF{time} + K{time}*(y(:,time) - z);
    Y_EKF_ALS_EKF{time} = yGT(:,time) - z;
    P{time} = P{time} - K{time}*C{time}*P{time};
end

EKF_ALS_EKF = xHatEKFPrior;

%MHE-ALS-EKF Error
[meanEKFErrors, maxEKFErrors, stdEKFErrors] = calcEKFError(EKF_ALS_EKF, xGT);
fprintf('Calculated errors for EKF-ALS-EKF filter: \n');
fprintf('mean: %e  \n', meanEKFErrors);
fprintf('max: %e \n', maxEKFErrors);
fprintf('std: %e \n', stdEKFErrors);

clearvars -except xGT yGT MHE_ALS_EKF Y_MHE_ALS_EKF Original_EKF Y_Original_EKF EKF_ALS_EKF Y_EKF_ALS_EKF

%% Convert to matrices
Original_EKF_mat(:,1) = Original_EKF{1,2};
EKF_ALS_EKF_mat(:,1) = EKF_ALS_EKF{1,2};
MHE_ALS_EKF_mat(:,1) = MHE_ALS_EKF{1,2};

for time = 2:size(Original_EKF,2)-1
    Original_EKF_mat(:,time) = Original_EKF{1,time+1};
    EKF_ALS_EKF_mat(:,time) = EKF_ALS_EKF{1,time+1};
    MHE_ALS_EKF_mat(:,time) = MHE_ALS_EKF{1,time+1};
end
xGT = xGT(:,2:end);

%% Comparison Results
for plotn = 1:2:12
    figure();
    hold on
    title(sprintf('State Predictions %d',plotn));
    xlabel('Time (sec)');
    
    if(plotn <= 6)
        ylabel('Position (m)');
    else
        ylabel('Angle (rad)');
    end
    
    xvals = 1:size(xGT,2);
    xvals = 0.1*xvals;
    plot(xvals, xGT(plotn,:), 'k--');
    plot(xvals, Original_EKF_mat(plotn,:), 'b');
    plot(xvals, EKF_ALS_EKF_mat(plotn,:), 'r');
    plot(xvals, MHE_ALS_EKF_mat(plotn,:), 'g');
    
    legend('True Value','Initial EKF','EKF-ALS-EKF','MHE-ALS-EKF');
    grid minor
    grid on
    title('Title here');
end
% tilefigs;

run plotCorrelations;

% %This figure is no longer used
% %Plot the frequency spectrum for the first two observations
% for n = 1:10
%     innovs_Orig = cell2mat(Y_Original_EKF);
%     innovs_EKF_ALS_EKF = cell2mat(Y_EKF_ALS_EKF);
%     innovs_MHE_ALS_EKF = cell2mat(Y_MHE_ALS_EKF);
%     
%     figure();
%     hold on;
%     Finnov1 = abs(fft(innovs_Orig(n,:)));
%     Finnov2 = abs(fft(innovs_EKF_ALS_EKF(n,:)));
%     Finnov3 = abs(fft(innovs_MHE_ALS_EKF(n,:)));
%     
%     %Normalise it
%     maxInnov = max([Finnov1(:);Finnov2(:);Finnov3(:)]);
%     Finnov1 = Finnov1./maxInnov;
%     Finnov2 = Finnov2./maxInnov;
%     Finnov3 = Finnov3./maxInnov;
%     
%     plot(Finnov1(1:floor(end/2)),'b');
%     plot(Finnov2(1:floor(end/2)),'r');
%     plot(Finnov3(1:floor(end/2)),'g');   
%     
% %     vline(20,'k--'); %This is the N value used in ALS
%     ylabel('Innovation Magnitude (Normalized)');
%     xlabel('Timesteps (100ms each)');
%     
%     legend('Initial EKF','EKF-ALS-EKF','MHE-ALS-EKF', 'ALS Window Length');
% end
