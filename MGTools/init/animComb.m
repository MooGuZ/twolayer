function anim = animComb(A)
% ANIMCOMB would combine animations with background value BG

% check input arguments
assert(numel(size(A)) == 3, ...
    'animation matrix should has and only has three dimensions.');

anim = min(A,[],3);

end