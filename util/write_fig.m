function write_fig(num,name)
figure(num)
fname = [name '.eps'];
print('-depsc',fname)
assert(~unix(['epstopdf ' fname]));
fprintf('Wrote %s.{eps,pdf}\n',name)