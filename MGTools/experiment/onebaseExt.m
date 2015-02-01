function [anim,objValue,b,z,normMean] = onebaseExt(videoPath,z)
% ONEBASE learn optimal complex base function according to one video sequence

% Program Setting
nEpoch = 17;
nAdapt = 700;
nInfer = 700;
stepInit  = 1e-1;
stepLimit = 1e-7;

% Read GIF file into video sequence
v = gif2anim(videoPath);
% Initialize Complex Base Function
b = complex(randn(size(v,1),1),randn(size(v,1),1));

% Normalize video data
normMean = mean(v,2);
v = bsxfun(@minus,v,normMean);

% Initialize Theta by random number
% Experiment shows the random initilization can avoid local minimum
if ~exist('z','var')
    beta  = randn(1,size(v,2));
    theta = pi * randn(1,size(v,2));
else
    beta  = abs(z);
    theta = angle(z);
end

% Polarize Complex Base Function
alpha = abs(b);
phi = angle(b);

% E-M Algo
error = calError(v,alpha,phi,beta,theta);
objValue = error(:)' * error(:);
for epoch = 1 : nEpoch
    % Adapting Complex Base Function
    step = stepInit;
    for i = 1 : nAdapt
        newAlpha = alpha - step * dAlpha(error,alpha,phi,beta,theta);
        newPhi = phi - step * dPhi(error,alpha,phi,beta,theta);
        newError = calError(v,newAlpha,newPhi,beta,theta);
        newObjValue = newError(:)' * newError(:);
        while(newObjValue > objValue)
            step = step / 2;
            if step < stepLimit, break; end                
            newAlpha = alpha - step * dAlpha(error,alpha,phi,beta,theta);
            newPhi = phi - step * dPhi(error,alpha,phi,beta,theta);
            newError = calError(v,newAlpha,newPhi,beta,theta);
            newObjValue = newError(:)' * newError(:);
        end
        if step < stepLimit, break; end
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
        newTheta = theta - step * dTheta(error,alpha,phi,beta,theta);
        newBeta = beta - step * dBeta(error,alpha,phi,beta,theta);
        newError = calError(v,alpha,phi,newBeta,newTheta);
        newObjValue = newError(:)' * newError(:);
        while(newObjValue > objValue)
            step = step / 2;
            if step < stepLimit, break; end
            newTheta = theta - step * dTheta(error,alpha,phi,beta,theta);
            newBeta = beta - step * dBeta(error,alpha,phi,beta,theta);
            newError = calError(v,alpha,phi,newBeta,newTheta);
            newObjValue = newError(:)' * newError(:);
        end
        if step < stepLimit, break; end
        beta = newBeta;
        theta = newTheta;
        error = newError;
        objValue = newObjValue;
    end
    disp(['Objective Value after infering process of EPOCH[',num2str(epoch), ...
        '] >> ',num2str(objValue)]);
end

b = alpha.*exp(1i*phi);
z = beta.*exp(1i*theta);

anim = bsxfun(@plus,normMean,real(repmat(b,1,size(v,2)) ...
    .* repmat(conj(z),size(v,1),1)));
% anim = real(repmat(b,1,size(v,2)) .* repmat(conj(z),size(v,1),1));
end

function e = calError(v,alpha,phi,beta,theta)
Theta = repmat(theta,size(v,1),1);
Beta  = repmat(beta,size(v,1),1);
Alpha = repmat(alpha,1,size(v,2));
Phi   = repmat(phi,1,size(v,2));

e = v - Alpha .* Beta .* cos(Phi - Theta);
end

function d = dAlpha(error,~,phi,beta,theta)
Theta = repmat(theta,length(phi),1);
Beta  = repmat(beta,length(phi),1);
Phi   = repmat(phi,1,length(theta));

d = -2 * sum(error .* Beta .* cos(Phi - Theta),2);
end

function d = dPhi(error,alpha,phi,beta,theta)
Theta = repmat(theta,length(phi),1);
Beta  = repmat(beta,length(phi),1);
Phi   = repmat(phi,1,length(theta));

d = 2 * alpha .* sum(error .* Beta .* sin(Phi - Theta),2);
end

function d = dBeta(error,alpha,phi,~,theta)
Theta = repmat(theta,length(phi),1);
Alpha = repmat(alpha,1,length(theta));
Phi   = repmat(phi,1,length(theta));

d = -2 * sum(error .* Alpha .* cos(Phi - Theta),1);
end

function d = dTheta(error,alpha,phi,beta,theta)
Theta = repmat(theta,length(phi),1);
Alpha = repmat(alpha,1,length(theta));
Phi   = repmat(phi,1,length(theta));

d = -2 * beta .* sum(error .* Alpha .* sin(Phi - Theta),1);
end
