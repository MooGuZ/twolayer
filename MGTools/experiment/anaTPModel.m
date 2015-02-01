function anaTPModel(M,output)
r = 32;
c = 32;

if ~exist(output,'dir')
    system(['mkdir -p ',output]);
end

f = figure('visible','off');

[npixel,npattern] = size(M.patBase);
disp('Draw Pattern Functions...');
for i = 1 : npattern
    I = mat2img(reshape(M.patBase(:,i),[r,c]), ...
        atan(tan(.45*pi)/max(abs(M.patBase(:,i))))/(pi/2));
    imwrite(I,[output,'/PatternBase-',num2str(i),'.png'],'PNG');
end

ntrans = size(M.transBase,2);
disp('Draw Transformation Functions...');
for i = 1 : ntrans
    I = mat2img(exp(1j*reshape(M.transBase(:,i),[r,c])));
    imwrite(I,[output,'/TransformBase-',num2str(i),'.png'],'PNG');
end

% disp('Animate Complex Base Functions...');
% B = complex(zeros(npixel,npattern*ntrans),zeros(npixel,npattern*ntrans));
% for i = 1 : npattern
%     for j = 1 : ntrans
%         B(:,(j-1)*npattern+i) = ...
%             M.patBase(:,i) .* exp(1j*M.transBase(:,j));
%     end
% end
% animateBase(B,'comb','gif','Filename',[output,'/ComplexBase'],'ncolor',256);

disp('Save Recovered Animation...');
ffindex = [M.ffindex,size(M.recoveredAnim,2)+1];
for i = 1 : length(ffindex)-1
    anim2gif(M.recoveredAnim(:,ffindex(i):ffindex(i+1)-1), ...
        [r,c,ffindex(i+1)-ffindex(i)],...
        [output,'/RecoveredAnim-',num2str(i),'.gif']);
end

close(f);
end