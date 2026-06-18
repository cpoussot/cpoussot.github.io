function saveGIF(handler,out,name,speed)

    if nargin < 4
        speed = .1;
    end
    F         = getframe(handler); 
    [RGB,~]   = frame2im(F); 
    [IND,map] = rgb2ind(RGB,255); 
    if out == 1
        imwrite(IND,map,[name '.gif'],'gif','LoopCount',Inf); 
    else
        imwrite(IND,map,[name '.gif'],'gif','WriteMode','append','DelayTime',speed); 
    end
end