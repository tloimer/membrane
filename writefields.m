fid = fopen('tmpstates','w');
for i = 1:length(psat1)
if ~isempty(fl0{i}), state0 = fl0{i}.sol.states; else state0 = []; end
if ~isempty(fl90{i}), state90 = fl90{i}.sol.states; else state90 = []; end
fprintf(fid,'%4u\t%s\t%s\n',i,state0,state90);
end
fclose(fid)
