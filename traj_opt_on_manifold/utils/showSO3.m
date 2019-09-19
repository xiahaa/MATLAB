function showSO3(Rdata,Rreg)
    Rdata = reshape(Rdata,3,3,[]);
    Rreg = reshape(Rreg,3,3,[]);
    N1 = size(Rdata,3);
    N2 = size(Rreg,3);
    xid = zeros(3,N1);
    xi = zeros(3,N2);
    for i = 1:N1
        xid(:,i)=logSO3(Rdata(:,:,i));
    end
    for i = 1:N2
        xi(:,i) = logSO3(Rreg(:,:,i));
    end
    v = [1;1;1]/sqrt(3);
    x1 = [xid(1,:).*v(1);xid(2,:).*v(2);xid(3,:)*v(3)];
    x2 = [xi(1,:).*v(1);xi(2,:).*v(2);xi(3,:)*v(3)];
    
    [x,y,z]=sphere(100);
    figure
    surf(x,y,z,'FaceAlpha',0.1,'FaceColor',[0.1,0.1,0.1],'EdgeColor','none')
    shading interp
    hold on;
    plot3(x1(1,:),x1(2,:),x1(3,:),'-o','MarkerSize',10,'Color','k');
    plot3(x2(1,:),x2(2,:),x2(3,:),'--s','MarkerSize',10,'Color','b');
end