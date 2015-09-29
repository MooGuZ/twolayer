function [a, b] = simsqrt(num)
% SIMSQRT factorize given number to two number as close as posible
% - this is a simple method can handle major situation, it, however, is not
% - a complete solution.
%
% USAGE : [A, B] = SIMSQRT(NUM)
%
% MooGu Z. <hzhu@case.edu>
% Jul 10, 2015

% CHAGE LOG
% Jul 10, 2015 - Version 0.00 : initial commit

% calculate factor list
facts = factor(num);
% ensure there are even number of factors in list
if mod(numel(facts), 2)
    facts = [1, facts];
end
% calculate the start position
spos = floor(numel(facts) / 4) + 1;
% calculate the end position
epos = spos + numel(facts) / 2 - 1;
% calculate <a> and <b>
a = prod(facts(spos : epos));
b = num / a;

end