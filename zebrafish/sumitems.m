function sumf = sumitems(force,pt,plist)
% sum up all items
sumf = 0;
if ~isempty(force) % normr normalize each row
    sumf = sum(repmat(force,1,2).*normr(plist-repmat(pt,size(plist,1),1))); % forces are row vectors, so we do the transpose 
end