function rec = objRecComb(recA,recB)
% This helper function combine two objective records into one record

% fundamental setting
prate  = 0.8; % rate of priority for later records
maxlen = 70; % maximal length of output records
% check length of each records
lenA = numel(recA.n);
lenB = numel(recB.n);
% cumulate filed 'n' of record B
for i = 2 : lenB, recB.n(i) = recB.n(i-1) + recB.n(i); end
% iteration count of each records
countA = recA.n(end);
countB = recB.n(end);
% combine two records
if lenA + lenB <= maxlen
    % shortcut for short records
    rec.n = [recA.n,recB.n+countA];
    rec.v = [recA.v,recB.v];
else
    % initialize output record
    rec.n = zeros(1,maxlen);
    rec.v = repmat(struct('noise',0,'sparse',0,'slow',0, ...
        'smpat',0,'smtrans',0,'value',0),1,maxlen);  
    % chose random sample according to iteration counts
    indA = unique(round(linspace(1,lenA, ...
        ceil((prate*maxlen*countA)/(countA+countB)))));
    indA = indA(1:end-1);
    ocpA = numel(indA);
    rec.n(1:ocpA) = recA.n(indA);
    rec.v(1:ocpA) = recA.v(indA);
    indB = unique(round(linspace(1,lenB,maxlen-ocpA)));
    ocpB = numel(indB);
    rec.n(ocpA+1:ocpA+ocpB) = recB.n(indB) + countA;
    rec.v(ocpA+1:ocpA+ocpB) = recB.v(indB);
    if ocpA + ocpB < maxlen
        rec.n = rec.n(1:ocpA+ocpB);
        rec.v = rec.v(1:ocpA+ocpB);
    end
end    

end