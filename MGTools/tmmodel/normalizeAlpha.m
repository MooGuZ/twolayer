function [alpha,phi,beta,theta,bia,delta,obj] = normalizeAlpha( ...
    alpha,phi,beta,theta,bia,delta,obj,sigma,ctrl,v,ffindex,res)

% initialization flag of recalculating objective value
recalc = false;

% define tolerance
tol = 1e-1;

if ctrl.swANorm
    if ~ctrl.swNegCut
        % normalize alpha
        normFactor = sum(abs(alpha),1) / ctrl.anorm;
        alpha  = bsxfun(@rdivide,alpha,normFactor);
        bia    = bsxfun(@times,bia,normFactor);
        beta   = bsxfun(@times,beta,normFactor);
        recalc = true;
    end
else
    % Scale factors of alpha
    a = max(abs(alpha));
    %Normalization and Transformation by Scale Factor
    if any((a-1) > tol)
        alpha = bsxfun(@rdivide,alpha,a);
        bia   = bsxfun(@times,bia,a);
        beta  = bsxfun(@times,beta,a);
        recalc = true;
    end
end

% Positive map of alpha
P = (alpha < 0);
% Set X and I
setX = any(P,2);
setI = any(P,1);
% Normalization and Transformation by Positive Restraint
if any(setX) || any(setI)
    alpha = abs(alpha);
    % Approximate Compensation of Absolute Operator
    if sum(setX)*sum(setI) > sum(~setX)*sum(~setI)
        phi(setX(:),:,:,:) = phi(setX(:),:,:,:) + pi;
        theta(:,~setI(:),:,:) = theta(:,~setI(:),:,:) + pi;
    end
%     % Send warning for disability of compansation normalization
%     if (sum(~setX) * sum(~setI) ~= 0) || ...
%             any(all([permute(any(abs(bia)>tol,4),[2,1,4,3]),setI(:)],2))
%         disp('Negative Flipping Method Applied!');
%     end
    % set objective value recalculation flag
    recalc = true;
end

% Recalculate Objective Value
if recalc
    delta = v - genmodel(alpha,phi,beta,theta,bia);
    obj   = objFunc(alpha,phi,beta,theta,bia,delta,sigma,ffindex,res);
    disp(['Objective Value after normalization of alpha >> ', ...
        num2str(obj.value)]);
end

end
