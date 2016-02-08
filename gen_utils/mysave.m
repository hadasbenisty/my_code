function mysave(f,filename)

saveas(f,[filename '.fig'],'fig');
saveas(f,[filename '.jpg'],'jpg');
print( f, '-depsc', [filename '.eps']);

