function drawBox(box,origin)
    L = box;
    xy = [ 
        0      0      0
        0      L(2)   0
        L(1)   L(2)   0
        L(1)   0      0
        0      0      0];
    
    xy = bsxfun(@plus, xy, origin.*box);
    
    hstate = ishold();
    plot3(xy(:,1),xy(:,2),xy(:,3),'-k'), hold on
    plot3(xy(:,1),xy(:,2),xy(:,3)+L(3),'-k')
    for i=1:5
        plot3(xy([i i],1),xy([i i],2),xy([i i],3)+[0 L(3)]','-k')
    end
    axis equal
    if ~hstate
        hold off
    end