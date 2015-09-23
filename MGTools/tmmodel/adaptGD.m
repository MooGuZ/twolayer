function [alpha,phi,delta,obj,i,ctrl] = ...
    adaptGD(niter,alpha,phi,beta,theta,bia,delta,obj, ...
        sigma,ctrl,v,ffindex,resolution)

% deal with easy case
if niter < 1 || ~(ctrl.swPatOpt || ctrl.swTransOpt)
    i = 0; return
end

% randomlize switcher of negative cutter
ctrl.swNegCut = rand(1) < ctrl.probNegCut;

step = ctrl.adaptInitStep;
for i = 1 : niter
    % calculate derivatives of bases
    [dAlpha,dPhi] = dBase(alpha,phi,beta,theta,bia, ...
        delta,sigma,ctrl,resolution);
    % step size search
    while(true)
        % update alpha
        if ctrl.swPatOpt
            newAlpha = alpha - step * dAlpha;
            % keep alpha positive
            if ctrl.swNegCut,
                newAlpha(newAlpha < 0) = 0;
                % normalize alpha
                if ctrl.swANorm
                    newAlpha = bsxfun(@rdivide,newAlpha, ...
                        sum(abs(newAlpha),1)/ctrl.anorm);
                end
            end
        else
            newAlpha = alpha;
        end
        % update phi
        if ctrl.swTransOpt
            newPhi = phi - step * dPhi;
        else
            newPhi = phi;
        end
        % calculate new objective value
        newDelta = v - genmodel(newAlpha,newPhi,beta,theta,bia);
        newObj   = ...
            objFunc(newAlpha,newPhi,beta,theta,bia,newDelta, ...
                sigma,ffindex,resolution);
        % update step size when test step does not 
        % - improve objective value
        if newObj.value > obj.value
            step = step / 2;
            if step < ctrl.accuracy
                return
            else
                continue
            end
        end
        % make the step
        alpha = newAlpha;
        phi   = newPhi;
        delta = newDelta;
        obj   = newObj;
        % go to next iteration
        break
    end
end

end
