function dispFirstlayer(m)

if isfield(m,'dewhitenMatrix')
    A = m.dewhitenMatrix*m.A;
    % dispBase(A,'dewhitenBase');
    A = m.zerophaseMatrix*A;
    dispBase(A,'zerophaseBase');
else
    dispBase(m.A,'orginBase');
end

end
