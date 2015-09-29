function arr = randpermext(n, k)
% RANDPERMEXT is the extention of <randperm> that generate random permutation of numbers
%
% [USAGE] arr = randpermext(n,k) generate random permutation of (1 to n), however, when k
%           is greater than n, then it generate a sequence of numbers that let each number
%           in (1 to n) appears almost same times with a random order.
%
% MooGu Z. <hzhu@case.edu>
% June 12, 2015 - Version 0.00 : initial commit

if exist('k', 'var')
    if k > n
        arr = [repmat(1:n, 1, floor(k/n)), randperm(n, rem(k,n))];
        arr = arr(randperm(k));
    else
        arr = randperm(n, k);
    end
else
    arr = randperm(n);
end

end