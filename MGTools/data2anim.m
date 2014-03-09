function frames = data2anim(m,p)
% DATA2ANIM read animation in column frames from data chunks specified by
% parameters

frames = crop_chunk(load_datachunk(m,p),m,p);

frames = bsxfun(@plus,m.dewhitenMatrix*frames,m.imageMean);

end