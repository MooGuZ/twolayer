function [anim,objValue,b,theta,normMean] = onebase(videoPath,theta)
% ONEBASE learn optimal complex base function according to one video sequence

% Program Setting
nEpoch = 13;
nAdapt = 100;
nInfer = 100;
stepInit  = 1e-1;
stepLimit = 1e-7;

% Read GIF file into video sequence
v = gif2anim(videoPath);
% Initialize Complex Base Function
b = complex(randn(size(v,1),1),randn(size(v,1),1));

% Normalize video data
normMean = mean(v,2);
v = bsxfun(@minus,v,normMean);

% Initilize Theta by average frequency estimation
% F = zeros(1,size(v,2));
% for i = 1 : size(v,1)
%     F = F + abs(fft(v(i,:)));
% end
% F = F / sum(F);
% if mod(size(v,2),2) == 0
%     theta = 2 * pi * ([0:(size(v,2)/2),-(size(v,2)/2)+1:-1] * F');
% else
%     theta = 2 * pi * ([0:(size(v,2)-1)/2,-(size(v,2)-1)/2:-1] * F');
% end
% Initialize Theta by random number
% Experiment shows the random initilization can avoid local minimum
if ~exist('theta','var')
    theta = pi * randn(1);
end

% Polarize Complex Base Function
alpha = abs(b);
phi = angle(b);

% E-M Algo
error = calError(v,alpha,phi,theta);
objValue = error(:)' * error(:);
for epoch = 1 : nEpoch
    % Adapting Complex Base Function
    step = stepInit;
    for i = 1 : nAdapt
        newAlpha = alpha - step * dAlpha(error,alpha,phi,theta);
        newPhi = phi - step * dPhi(error,alpha,phi,theta);
        newError = calError(v,newAlpha,newPhi,theta);
        newObjValue = newError(:)' * newError(:);
        while(newObjValue > objValue)
            step = step / 2;
            if step < stepLimit break; end                
            newAlpha = alpha - step * dAlpha(error,alpha,phi,theta);
            newPhi = phi - step * dPhi(error,alpha,phi,theta);
            newError = calError(v,newAlpha,newPhi,theta);
            newObjValue = newError(:)' * newError(:);
        end
        if step < stepLimit break; end
        alpha = newAlpha;
        phi = newPhi;
        error = newError;
        objValue = newObjValue;
    end
    disp(['Objective Value after adapting process of EPOCH[',num2str(epoch), ...
        '] >> ',num2str(objValue)]);
    % Infering Optimal Theta
    step = stepInit;
    for i = 1 : nInfer
        newTheta = theta - step * dTheta(error,alpha,phi,theta);
        newError = calError(v,alpha,phi,newTheta);
        newObjValue = newError(:)' * newError(:);
        while(newObjValue > objValue)
            step = step / 2;
            if step < stepLimit 
                break; 
            end
            newTheta = theta - step * dTheta(error,alpha,phi,theta);
            newError = calError(v,alpha,phi,newTheta);
            newObjValue = newError(:)' * newError(:);
        end
        if step < stepLimit 
            break; 
        end
        theta = newTheta;
        error = newError;
        objValue = newObjValue;
    end
    disp(['Objective Value after infering process of EPOCH[',num2str(epoch), ...
        '] >> ',num2str(objValue)]);
end

b = alpha.*exp(1i*phi);

anim = bsxfun(@plus,normMean,real(repmat(b,1,size(v,2)) ...
    .* repmat(exp(-1i*theta*(0:size(v,2)-1)),size(v,1),1)));

end

function e = calError(v,alpha,phi,theta)
T = repmat(0:size(v,2)-1,size(v,1),1);
Alpha = repmat(alpha,1,size(v,2));
Phi = repmat(phi,1,size(v,2));

e = v - Alpha .* cos(Phi - theta * T);
end

function d = dAlpha(error,~,phi,theta)
T = repmat(0:size(error,2)-1,size(error,1),1);
Phi = repmat(phi,1,size(error,2));

d = -2 * sum(error .* cos(Phi - theta * T),2);
end

function d = dPhi(error,alpha,phi,theta)
T = repmat(0:size(error,2)-1,size(error,1),1);
Alpha = repmat(alpha,1,size(error,2));
Phi = repmat(phi,1,size(error,2));

d = 2 * sum(error .* Alpha .* sin(Phi - theta * T),2);
end

function d = dTheta(error,alpha,phi,theta)
T = repmat(0:size(error,2)-1,size(error,1),1);
Alpha = repmat(alpha,1,size(error,2));
Phi = repmat(phi,1,size(error,2));

d = -2 * sum(sum(error .* Alpha .* sin(Phi - theta * T) .* T));
end
