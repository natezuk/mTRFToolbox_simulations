function clridx = get_color_idx(idx,max_idx,cmap)
% compute the index in cmap that allows an even distribution of colors
% between 0 and max_idx (for example, from a for loop) along cmap

clridx = round((idx-1)/max_idx*(size(cmap,1)-1))+1;