function funPlotGLRAfreqWing(w,Hf,H,name)
    
    r       = length(H.a);
    nw      = 2e3;
    w_thin  = union(w,logspace(log10(min(w)),log10(max(w)),nw));
    el      = 30;
    az      = 120;
    Hdata   = mdspack.freqresp(Hf,w);
    Hmdl    = mdspack.freqresp(H,w_thin);
    
    handler=figure('Color','white'); hold on;
    col         = get(gca,'colororder');
    %subSetOut   = ([23:2:31 33:44]);
    subSetOut   = 1:size(H,1);
    subSetImage = flip([1:12 14:2:22]);
    nPoints     = numel(subSetOut);
    location    = linspace(0,25,nPoints);%-offset/2-.5;
    for out = 1:nPoints 
        subplot(3,3,[1 2 4 5 7 8]), hold on, grid on
        plot3(location(out)*ones(numel(w'),1),w',mag2db(abs(squeeze(Hdata(subSetOut(out),1,:)))),'.','Color',col(1,:),'MarkerSize',25)
        Z(:,out) = mag2db(abs(squeeze(Hmdl(subSetOut(out),1,:))));
        Y(:,out) = w_thin';
        X(:,out) = location(out)*ones(numel(w_thin'),1);
        %axis tight
        zlim([-65 -20])
        xlabel('Wingspan','interpreter','latex','FontSize',20); 
        ylabel('Frequency [rad/s]','interpreter','latex','FontSize',20); 
        zlabel('Gain [dB]','interpreter','latex','FontSize',20)
        set(gca,'YScale','log')
        names = {'Wing root';' ';'Wing tip'};
        xlim([1 location(nPoints)])
        set(gca,'xtick',[1 17 24],'xticklabel',names,'TickLabelInterpreter','latex');
        set(gca,'ytick',[1 10 50 100],'xticklabel',names,'TickLabelInterpreter','latex');
    end
    surfc(X,Y,Z,'FaceColor',[1 1 1]*.9,'EdgeColor','none','DisplayName','Identified model'), 
    title(['Frequency response surface along wingspan ($r=' num2str(r) '$)'],'interpreter','latex','FontSize',20)
    
    kk = 0;
    for out = 1:nPoints 
        subplot(3,3,[1 2 4 5 7 8]), cla,
        surf(X,Y,Z,'FaceColor',[1 1 1]*.9,'EdgeColor','none'),
        if out < nPoints
            surf(X(:,out:(out+1)),Y(:,out:(out+1)),Z(:,out:(out+1)),'EdgeColor','none','DisplayName','Identified model')
        else
            surf(X,Y,Z,'EdgeColor','none','DisplayName','Identified model'),
        end
        for slk = 1:nPoints
            plot3(location(slk)*ones(numel(w'),1),w',mag2db(abs(squeeze(Hdata(subSetOut(slk),1,:)))),'.','Color',col(1,:),'MarkerSize',25,'DisplayName','CFD data')
        end
        %legend('show','location','northeast')
        view(az+out*2,el)
        %
        subplot(3,3,[3 6])
        img = imread(['CutJPG/Cut-' num2str(subSetImage(out),'%02i') '.jpg']);
        image(img);
        set(gca,'XLim',[200 1500])
        set(gca,'YLim',[800 1600])
        set(gca,'xtick',[],'ytick',[])
        title('($10^5$ nodes, Mach 0.84, $U_\infty=257.93$m/s)','interpreter','latex','FontSize',20)
        %set(gca,'visible','off');
        %
        % subplot(2,3,6)
        % img = imread('figures/CpA2.png');
        % image(img);
        % set(gca,'xtick',[],'ytick',[])
        % set(gca,'visible','off');
        %
        drawnow
        kk = kk + 1;
        saveGIF(handler,kk,name)
        %mdspack.plot.figSavePDF(['flexibleAC_' num2str(kk)],.5)
    end
    
    rot0    = az+out*2;
    rot     = linspace(rot0,360+az,50);
    subplot(3,3,[3 6]), cla
    img = imread('CpA2.jpeg');
    image(img);
    set(gca,'xtick',[],'ytick',[])
    %set(gca,'visible','off');
    title('($10^5$ nodes, Mach 0.84, $U_\infty=257.93$m/s)','interpreter','latex','FontSize',20)
    subplot(3,3,[1 2 4 5 7 8]), cla,
    for slk = 1:nPoints
        plot3(location(slk)*ones(numel(w'),1),w',mag2db(abs(squeeze(Hdata(subSetOut(slk),1,:)))),'.','Color',col(1,:),'MarkerSize',25)
    end
    surf(X,Y,Z,'EdgeColor','none'),
    alpha(.7)
    for ii = rot
        view(ii,el)
        kk = kk + 1;
        drawnow
        saveGIF(handler,kk,name)
    end

end
%%%
function saveGIF(handler,out,name)
    F         = getframe(handler); 
    [RGB,~]   = frame2im(F); 
    [IND,map] = rgb2ind(RGB,255); 
    if out == 1
        imwrite(IND,map,[name '.gif'],'gif','LoopCount',Inf); 
    else
        imwrite(IND,map,[name '.gif'],'gif','WriteMode','append','DelayTime',.1); 
    end
end