function clef()
% closes empty figures
fig_array = get(0, 'Children');
for i = 1 : numel(fig_array)
   if isempty(get(fig_array(i), 'Children'))
       close(fig_array(i));
   end
end