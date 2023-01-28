

mina = 2000;
maxa = 7000;
da = 500;

minb = 300;
maxb = 500;
db = 10;

for a=mina:da:maxa
    for b=minb:db:maxb
        filename = ['data/run2_a',num2str(a),'_b',num2str(b),'.mat'];
        load(filename)
        plt = size(S,2);
        ntype =2;
        rall=[20,20];
        colorlist = {[0 0 0],[1 0.65 0],[1 1 0],[0.65 0.67 1],[0.67 0.78 0.94]}; % black, orange, yellow, light blue, white
        translist = [0.7 0.3 0.3 0.3 0.3]; % list of transparency value
        for indt = plt
            %hold on
        %     f1 = axes('Position',[0.1 0.1 0.77 0.81]);
        %     hold(f1,'on');
            hold on
            for indi = 1:ntype
                plotcells(S(indt).pos{indi},rall(indi),colorlist{indi},translist(indi));
                title(['Time at ', num2str(indt),' days']);
                plot([0 S(indt).domsize(1)],[0 0],'b');
                plot([0 S(indt).domsize(1)],[S(indt).domsize(2) S(indt).domsize(2)],'b');
                axis equal
                %axis([0 5000 0 4000]);
            end
            % hold(f1,'off');
%             if indt ~= plt
%                 clf;
%             end
            figname =  ['figures/figs_d51_a',num2str(a),'_b',num2str(b),'.png'];
            saveas(gcf,figname);
            clf
            clearvars S
        end

    end
end