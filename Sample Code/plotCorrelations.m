%Based on innovation data, plot the correlation of the lags over time

clearvars
% close all
load('innovations.mat');

figure()
hold on

%For each observation element
for n = 1:5
    dataMat = [];
    maxLags = 90;

    %Need to have a matrix where each column is N lags back, each row is
    %observations from that number of lags
    for k = 10:500
        newRow = innovs_Orig(n,k:k+maxLags);
       dataMat = [dataMat; newRow]; 
    end

    R = corrcoef(dataMat);

    plot(1:size(R,1), R(1,:));
    title('Innovations Correlation Over Time');
    xlabel('Time Steps');
    ylabel('Correlation Coefficient');
    grid minor
    grid on
    
end
tilefigs

