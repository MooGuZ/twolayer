function dispMCode(C,crds,tstr)

% Parameters Setting
clr   = [.5,.5,.5]; % Color of Motion Vector Arrow
alpha = 0.07;       % Scale Factor of Motion Vector Arrow
gamma = 2.4;        % Gamma Correction

% Gamma Correction Function
gcorrect = @(I) I.^gamma;

swPrint = exist('tstr','var') && ~isempty(tstr);

[nbase,ncode] = size(C);

% Initialize Figure Window
if swPrint
    figure('visible','off');
else
    figure();
end
set(gcf, ...
    'Color',clr, ...
    'inverthardcopy','off', ...
    'position',[100,100,1024,1024]);

% Calculate Spacial Frequency of First Layer Bases
spfreq = sqrt(sum(crds(3:4,:).^2));
% Check Constant Component and Remove Corresponding Code Component
if any(spfreq==0)
    C(spfreq==0,:) = 0;
    spfreq = spfreq + eps;
end
% Motion Vector Direction for Each Bases
mvdir = bsxfun(@rdivide,crds(3:4,:),spfreq);

for j = 1 : ncode
	% Calculate Motion Vector
    mvamp = alpha * C(:,j)' ./ spfreq;
    xmvec = mvamp .* mvdir(1,:);
    ymvec = mvamp .* mvdir(2,:);
    % Sort according to Phase Difference
    [~,index] = sort(abs(C(:,j)),'ascend');
    % Calculate Colormap according to Phase Change
    ineg = C(index,j) < 0;
    cmap = gcorrect(abs(C(index,j))/pi);
    cmap(ineg,j) = -cmap(ineg,j);
    cmap = repmat((cmap+1)/2,1,3);
    % Draw Center Points of Bases
    if swPrint
        scatter(crds(1,index),crds(2,index),17,cmap,'filled');
    else
        scatter(crds(1,index),crds(2,index),70,cmap,'filled');
    end
    axis equal off
    hold on
    % Draw Motion Vector
    if swPrint
        for i = 1 : nbase
            quiver(crds(1,index(i)),crds(2,index(i)), ...
                xmvec(index(i)),ymvec(index(i)),0, ...
                'Color',cmap(i,:));
        end
    else
        for i = 1 : nbase
            quiver(crds(1,index(i)),crds(2,index(i)), ...
                xmvec(index(i)),ymvec(index(i)),0, ...
                'Color',cmap(i,:),'LineWidth',1.5);
        end
    end
    colormap gray
    colorbar
    % Figure Setting
    hold off
    % Print
    if swPrint
        if ncode > 1
            print(gcf,'-depsc2','-r300',...
                [tstr,'[',num2str(j),'].eps']);
        else
            print(gcf,'-depsc2','-r300',...
                [tstr,'.eps']);
        end
    end
end

end
