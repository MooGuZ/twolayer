function strarr = strsplitby(str, delimiter)
% STRSPLITBY split string by specified delimiter
%
% [USAGE] strarr = strsplitby(str, delimiter)
%
% MooGu Z. <hzhu@case.edu>
% June 12, 2015 - Version 0.00 : initial commit

assert(ischar(str), 'First argument needs to be a string.');
assert(ischar(delimiter) && numel(delimiter) == 1, ...
    'Delimiter (second argument) needs to be a character.');

% search for delimiter
index = [0, find(str == delimiter, inf), numel(str) + 1];
% create a cell array of string
strarr = cell(0);
for i = 1 : numel(index) - 1
    substr = str(index(i) + 1 : index(i+1) - 1);
    if ~isempty(substr)
        strarr = [strarr, {substr}];
    end
end

end