function [meanErrors, maxErrors, stdErrors] = calcEKFError(states, groundTruth)
    
clearvars delta
if(iscell(states))
    states = cell2mat(states(2:end));
end

deadzone = 50; %Ignore starting / ending points
delta = states(:,deadzone:(end-deadzone)) - groundTruth(:,deadzone+1:(end-deadzone));
deltaMax = states(:,deadzone:(end-deadzone)) - groundTruth(:,deadzone+1:(end-deadzone));

% delta = states-groundTruth(:,2:end);
% deltaMax = states(:,10:end) - groundTruth(:,11:end);

rn = 1;
for row = 1:2:12
    meanErrors(rn) = sqrt(mean(delta(row,:).^2)); 
    maxErrors(rn) = max(abs(deltaMax(row,:)));
    stdErrors(rn) = std(delta(row,:));
    rn = rn + 1;
end

end
