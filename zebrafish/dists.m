function ds = dists(npt, plist)
if ~isempty(plist) && ~isempty(npt)
    ds = sqrt((npt(1)-plist(:,1)).^2 + (npt(2)-plist(:,2)).^2); % column vector?
else
    ds = [];
end