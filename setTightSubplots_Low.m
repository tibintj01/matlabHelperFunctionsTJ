make_it_tight = true;
%subplot = @(m,n,p) subtightplot (m, n, p, [0.075 0.05], [0.05 0.05], [0.05 0.05]);
%function h=subtightplot(m,n,p,gap,marg_h,marg_w,varargin)
subplot = @(m,n,p) subtightplot (m, n, p, [0.1 0.1], [0.05 0.05], [0.1 0.1]);

if ~make_it_tight,  clear subplot;  end
