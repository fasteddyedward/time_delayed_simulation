function Hello_world(a,b)

f=figure;
f.Visible='off';
title([num2str(a),'Hello world!',num2str(b)])
saveas(gcf,['test_figure.png'])   
exit
end