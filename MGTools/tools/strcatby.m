function str = strcatby(strarr, symbol)
% STRCATBY concatenate strings in cell array with specified symbol
%
% [USAGE] str = strcatby(strarr, symbol)
%
% MooGu Z. <hzhu@case.edu>
% June 12, 2015 - Version 0.00 : initial commit

assert(iscell(strarr), 'First argument needs to be a cell array of string.');

str = strarr{1};
for i = 2 : numel(strarr)
    str = [str '.' strarr{i}];
end

end