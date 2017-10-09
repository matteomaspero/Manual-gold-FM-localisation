function printpdf(h,outfilename)
%printpdf Generate a pdf with high resolution
% This is ideal to produce pdf files suitable for publication!
dpi = 300;
set(h, 'PaperUnits','centimeters');
set(h, 'Units','centimeters');
pos=get(h,'Position');
set(h, 'PaperSize', [pos(3) pos(4)]);
set(h, 'PaperPositionMode', 'manual');
set(h, 'PaperPosition',[0 0 pos(3) pos(4)]);
print('-dpdf',outfilename,sprintf('-r%d',dpi));

end


