make_it_tight = true;
%function h=subtightplot(m,n,p,gap,marg_h,marg_w,varargin)
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.005], [0.1 0.1], [0.05 0.05]);
%subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.005], [0.1 0.1], [0.5 0.01]);
if ~make_it_tight,  clear subplot;  end
