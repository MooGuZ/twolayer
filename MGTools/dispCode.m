function dispCode(C,crds,tstr)

swPrint = exist('tstr','var');

[~,ncode] = size(C);

f = figure();
set(gcf,'Color','k');
set(gcf,'inverthardcopy','off');

for j = 1 : ncode
    figure(f); clf;
	[~,index] = sort(abs(C(:,j)),'descend');
	val = C(index,j) / max(abs(C(:,j)));
	% Spacial Domain
	subp(1,2,1);
	scatter(crds(1,index),crds(2,index),60,val,'filled');
	axis equal off;
	% title(['Code [',num2str(j),'] - Spacial Domain'],'Color','w');
	% Spacial Frequency Domain
	subp(1,2,2);
	scatter(crds(3,index),crds(4,index),60,val,'filled');
	axis equal off;
	% title(['Code [',num2str(j),'] - Spacial Freq Domain'],'Color','w');
    % Print
    if swPrint
        print(gcf,'-depsc2','-r300',...
            ['./fig/',tstr,'[',num2str(j),'].eps']);
    end
end

end
