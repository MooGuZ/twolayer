function m = randNormedMatrix(nrows, ncols)
% RANDNORMEDMATRIX works as a helper function to create random matrix with
% specified size that normalized in columns by 2-Norm equals to 1.
%
% USAGE : m = randNormedMatrix(nrows, ncols)
%
% MooGu Z. <hzhu@case.edu>
% June 16, 2015 - Version 0.00 : initial commit

m = randn(nrows, ncols);
m = bsxfun(@divide, m, sqrt(sum(m.^2)));

end