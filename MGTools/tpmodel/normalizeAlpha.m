function [alpha,phi,beta,theta,bia,delta,obj] = ...
    normalizeAlpha(alpha,phi,beta,theta,bia,delta,obj,res,v,sigma,ffindex)
zeroDecisionBound = 1e-2;
% Scale factors of alpha
a = max(abs(alpha));
% % Positive map of alpha
% P = (alpha < 0);
% % Set X and I
% setX = any(P,2);
% setI = any(P,1);
% Normalization and Transformation by Scale Factor
if any(abs(a-1) > zeroDecisionBound)
    alpha = bsxfun(@rdivide,alpha,a);
    bia   = bsxfun(@times,bia,a);
    beta  = bsxfun(@times,beta,a);
    % Recalculate Objective Value
    delta = v - genmodel(alpha,phi,beta,theta,bia);
    obj   = objFunc(alpha,phi,beta,theta,bia,delta,sigma,ffindex,res);
    disp(['Objective Value after normalization of alpha >> ', ...
        num2str(obj.value)]);
end
% % Normalization and Transformation by Positive Restraint
% if any(setX) || any(setI)
%     alpha = abs(alpha);
%     % Approximate Compensation of Absolute Operator
%     if sum(setX)*sum(setI) > sum(~setX)*sum(~setI)
%         phi(setX(:),:,:,:) = wrapToPi(phi(setX(:),:,:,:) + pi);
%         theta(:,:,~setI(:),:) = wrapToPi(theta(:,:,~setI(:),:) + pi);
%     end
%     % Send warning for disability of compansation normalization
%     if (sum(~setX) * sum(~setI) ~= 0) || ...
%             any(all([permute(any(abs(bia)>zeroDecisionBound,2),[3,1,2,4]), ...
%             setI(:)],2))
%         disp('[Warning] The effects of alpha normalization cannot be fully compansated!');
%     end
% end

end
