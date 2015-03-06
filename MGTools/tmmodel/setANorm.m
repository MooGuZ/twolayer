function model = setANorm(model,val)
% This is a temporal helper function to set switcher of alpha-normalization

% deal with special cases
if model.ctrl.swANorm == val, return
elseif val == false && isfield(model.ctrl,'anorm')
    model.ctrl = rmfiled(model.ctrl,'anorm');  
end
% set switcher to specified value
model.ctrl.swANorm = val;

end