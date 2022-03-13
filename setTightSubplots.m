make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.01], [0.05 0.01], [0.05 0.01]);
if ~make_it_tight,  clear subplot;  end
