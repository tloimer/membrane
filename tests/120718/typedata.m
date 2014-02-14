isdata = eval(identifikation);

% Drucke die Daten.
disp(sprintf(['#  ' identifikation '\n\n']));
kopf = '#\tkap/kl\tm/mgas\tm/mcalc\texp_id     \tP1\tp1-p2\tCcc\tT/Tcalc';
zeile = '%3u\t%5.3f\t%5.3f\t%5.3f\t%-11s\t%5.3f\t%4.1f\t%6.3f\t%5.3f';
disp(sprintf(kopf));
lendata = find(isdata);
for ii = 1:size(lendata,1)
  i = lendata(ii);
  disp(sprintf(zeile,i,Kk(i),meg(i),mec(i),...
    exp_id{i},P1(i),(p1(i)-p2(i))/1e3,Ccc(i),T12meas(i)/T12calc(i)));
end
