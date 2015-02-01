function dispBase(A,tstr)

swPrint = exist('tstr','var');

[npixel,nbase]=size(A);
sz=sqrt(npixel);

% Create Figure
if swPrint
    f = figure('visible','off');
else
    f = figure();
end
% Draw Bases one by one
for j = 1 : nbase
    sfigure(f); clf;
    % the j-th Base
    Base = reshape(A(:,j),sz,sz) / max(abs(A(:,j)));
    % Real Components
    subp(2,2,1);
    imagesc(real(Base),[-1,1]);
    colormap gray; freezeColors; axis image off;
    title(['Base[',num2str(j),'] - Real Components']);
    % Imaginary Components
    subp(2,2,2);
    imagesc(imag(Base),[-1,1]);
    colormap gray; freezeColors; axis image off;
    title(['Base[',num2str(j),'] - Imaginary Components']);
    % Pattern Components (Amplitude)
    subp(2,2,3);
    imagesc(abs(Base),[0,1]);
    colormap gray; freezeColors; axis image off;
    title(['Base[',num2str(j),'] - Pattern Components (Amplitude)']);
    % % Motion Components (Angle)
    % subp(2,2,4);
    % imagesc(angle(Base),[-pi,pi]);
    % colormap hsv; freezeColors; axis image off;
    % title(['Base[',num2str(j),'] - Motion Components (Angle)']);
    % Motion Components (Angle Modulated by Amplitude)
    subp(2,2,4);
    imshow(comp2img(Base));
    title(['Base[',num2str(j),'] - Motion Components (Phase<+Amp>)'])
    % Print Figure
    if swPrint
        if (nbase ~= 1)
            print(gcf,'-dpng',...
                [tstr,'[',num2str(j),'].png']);
        else
            print(gcf,'-dpng',[tstr,'.png']);
        end
    else
        f = figure();
    end
end

close(f);

end