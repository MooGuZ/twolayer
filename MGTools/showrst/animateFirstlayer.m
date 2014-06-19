function animateFirstlayer(m)
% ANIMATEFIRSTLAYER write animation of bases in first layer into GIF files

A = double(m.A);
% Dewhiten Bases if necessary
if isfield(m,'dewhitenMatrix')
    A = m.zerophaseMatrix*m.dewhitenMatrix*A;
end
% Write animation into GIF files
animateBase(A,'real','gif','Filename','zerophaseBase');

end