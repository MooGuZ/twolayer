function objplot(obj,swPrint,fname)

% assign default value
swPrint = exist('swPrint','var') && swPrint;
% check input arguments
assert(isstruct(obj), ...
    'input argument should be a structure of objective records!');
if swPrint
    assert(logical(exist('fname','var')),'File name is needed!');
end

% initialize a new figure
if swPrint
    f = figure('visible','off');
else
    f = figure();
end
% plot objective value
h(1) = semilogy(obj.n,[obj.v(:).value],'r.','markersize',13); hold on
semilogy(obj.n,[obj.v(:).value],'k-o','markersize',7); grid on
% plot sparseness
h(2) = semilogy(obj.n,[obj.v(:).sparse],'g.','markersize',13);
semilogy(obj.n,[obj.v(:).sparse],'k-o','markersize',7);
% plot slowness
h(3) = semilogy(obj.n,[obj.v(:).slow],'b.','markersize',13);
semilogy(obj.n,[obj.v(:).slow],'k-o','markersize',7);
% plot smoothness of pattern bases
h(4) = semilogy(obj.n,[obj.v(:).smpat],'y.','markersize',13);
semilogy(obj.n,[obj.v(:).smpat],'k-o','markersize',7);
% plot smoothness of transform bases
h(5) = semilogy(obj.n,[obj.v(:).smtrans],'c.','markersize',13);
semilogy(obj.n,[obj.v(:).smtrans],'k-o','markersize',7); hold off

% add label to axis
xlabel('#iteration');
ylabel('Objective Value');
title('Learning Curve');

% add legend
legend(h,'ObjValue','Sparseness','Slow-Changing', ...
    'PatSmoothness','TransSmoothness');

% print to figure
if swPrint
    print(f,'-depsc2','-r300',fname);
    % reset visibility and close figure
    set(f,'visible','on');
    close(f);
end

end

