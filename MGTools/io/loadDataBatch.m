function Data = loadDataBatch(m,p)
% LOADDATABATCH would load data into frame vectors according to setting in
% paramters
%
%   Data = loadDataBatch(m,p)
%
% MooGu Z. <hzhu@case.edu>
% June 17, 2014 - Version 0.1

switch p.data.type
    case 'fvp'
        Data = zeros(m.patch_sz^2,p.data.scope*p.data.nframe,'uint8');
        if p.data.scope >= p.data.quantity
            for i = 1 : p.data.scope
                Data(:,(i-1)*p.data.nframe+1:i*p.data.nframe) ...
                    = im2uint8(gif2anim( ...
                        [p.data.path,p.data.nameList{i}], ...
                        1:p.data.nframe));
            end
        else
            % Random Number Genereating
            s = randi(p.data.quantity,p.data.scope,1);
            % Load GIF Files
            for i = 1 : p.data.scope
                Data(:,(i-1)*p.data.nframe+1:i*p.data.nframe) ...
                    = im2uint8(gif2anim( ...
                        [p.data.path,p.data.nameList{s(i)}], ...
                        1:p.data.nframe));
            end
        end
end

end
        
