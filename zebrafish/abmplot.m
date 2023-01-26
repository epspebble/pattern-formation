function abmplot(S,ntype,totd)

figure()
plt = totd;
colorlist = {[0 0 0],[1 0.65 0],[1 1 0],[0.65 0.67 1],[0.67 0.78 0.94]}; % black, orange, yellow, light blue, white
translist = [0.7 0.3 0.3 0.3 0.3]; % list of transparency value

for indt = 1:plt
    hold on
    for indi = 1:ntype
        plotcells(S(indt).pos{indi},rall(indi),colorlist{indi},translist(indi));
        title(['Time at ', num2str(indt),' days']);
        plot([0 S(indt).domsize(1)],[0 0],'b');
        plot([0 S(indt).domsize(1)],[S(indt).domsize(2) S(indt).domsize(2)],'b');
        axis equal
        %axis([0 5000 0 4000]);
    end
    pause(0.2)
    hold off
    if indt ~= plt
        clf;
    end
end
end