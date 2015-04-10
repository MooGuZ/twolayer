function [alpha,phi,beta,theta,bia,delta,obj] = normalizeAlpha( ...
    alpha,phi,beta,theta,bia,delta,obj,sigma,ctrl,v,ffindex,res,verbose)

% initialization flag of recalculating objective value
recalc = false;

% define tolerance
tol = 1e-1;

if ctrl.swANorm
    if ~ctrl.swNegCut
        % calculate scale factors of alpha
        normFactor = sum(abs(alpha),1) / ctrl.anorm;
		% alpha-normalization and compansative transformation
        alpha  = bsxfun(@rdivide,alpha,normFactor);
        bia    = bsxfun(@times,bia,normFactor);
        beta   = bsxfun(@times,beta,normFactor);
		% set up recalculation flag
        recalc = true;
    end
else
    % calculate scale factors of alpha
    normFactor = max(abs(alpha));
    % alpha-normalization and compansative transformation
    if any((normFactor-1) > tol)
        alpha = bsxfun(@rdivide,alpha,normFactor );
        bia   = bsxfun(@times,bia,normFactor );
        beta  = bsxfun(@times,beta,normFactor );
		% set up recalculation flag
        recalc = true;
    end
end

% positive map of alpha
P = (alpha < 0);
% set X and I
setX = any(P,2);
setI = any(P,1);
% alpha-absolutization and compansative transformation
if any(setX) || any(setI)
    alpha = abs(alpha);
    % approximate compensation of absolute operator
    if sum(setX)*sum(setI) > sum(~setX)*sum(~setI)
        phi(setX(:),:,:,:) = phi(setX(:),:,:,:) + pi;
        theta(:,~setI(:),:,:) = theta(:,~setI(:),:,:) + pi;
    end
%     % Send warning for disability of compansation normalization
%     if (sum(~setX) * sum(~setI) ~= 0) || ...
%             any(all([permute(any(abs(bia)>tol,4),[2,1,4,3]),setI(:)],2))
%         disp('Negative Flipping Method Applied!');
%     end
    % set up recalculation flag
    recalc = true;
end

% Recalculate Objective Value
if recalc
    delta = v - genmodel(alpha,phi,beta,theta,bia);
    obj   = objFunc(alpha,phi,beta,theta,bia,delta,sigma,ffindex,res);
    if verbose >= 2
        disp(['Objective Value after normalization of alpha >> ', ...
            num2str(obj.value)]);
    end
end

end
