function qend = qend(zevap,m,f2ph,zend)
y2ph=deval(f2ph,zevap);
fg = intg(m,y2ph(1),ps(y2ph(1)),qevap(m,y2ph(1),y2ph(2)),[zevap zend]);
qend = fg.y(3,end);
