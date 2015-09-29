function value = findProperty(arglist, key, validfunc)
% FINDPROPERTY fetch value corresponds to the key in argument list
%
% VALUE = FINDPROPERTY(ARGLIST, KEY, VALIDFUNC)
%
% MooGu Z. <hzhu@case.edu>
% Jun  3, 2015 - Version 0.00 : initial commit
% Jun 25, 2015 - Version 0.01 : ignore all [space] in <key>

% flag of exist validation function
flagValidFunc = exist('validfunc', 'var');
if flagValidFunc    
    assert(isa(validfunc,'function_handle'), ...
        'the 3rd argument need to be a function handle for validation');
end

index = 1;
nargs = numel(arglist);
% traverse argument list to pick out keys
while index <= nargs
    arg = arglist{index};
    assert(ischar(arg), 'Property list is illegal.');
    % remove all [space] in <arg>
    arg = arg(arg ~= ' ');
    
    if arg(1) == '-' % unary (switcher) property
        if strcmpi(key, arg)
            value = true;
            return
        end
        % move to next one
        index = index + 1;
    else % binary (key value pair) property
        if strcmpi(key, arg)
            assert(nargs >= index + 1, ...
                sprintf('value of binary property [%s] is missing.', arg));
            value = arglist{index + 1};
            % validation check
            if flagValidFunc
                assert(validfunc(value), sprintf('Validation process failed for [%s]', key));
            end
            return
        end
        % move to next one
        index = index + 2;
    end
end

% don't match the key in argument list
if key(1) == '-' % key of unary property
    value = false;
else % key of binary property
    value = nan;
end
    
end
    