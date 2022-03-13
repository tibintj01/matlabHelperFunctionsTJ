h_fig = figure;
set(h_fig,'KeyPressFcn',@myfun);

function myfun(src,event)
   disp([event.Key  event.Key]);
end