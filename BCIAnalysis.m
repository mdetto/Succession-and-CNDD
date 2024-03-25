
clear
    gx1=[];gx2=[];
    gy1=[];gy2=[];
    dbh1=[];dbh2=[];
    ID1=[];ID2=[];
    c1=[];c2=[];
    agb1=[];agb2=[];
    yr = [1982 1985 1990 1995 2000 2005 2010 2015];  
for i=1:8
    
    
    dat = load(['C:\Users\mdetto\Dropbox (Princeton)\bci_Matlab_Files\bci' num2str(i) '.mat']);
    
    use1 = strcmp(dat.sp,'ceibpe');
    use2 = strcmp(dat.sp,'diptpa');
    
    gx1=[gx1;dat.gx(use1)];
    gy1=[gy1;dat.gy(use1)];
    dbh1=[dbh1;dat.dbh(use1)];
    agb1=[agb1;dat.agb(use1)];
    ID1=[ID1;dat.treeID(use1)];
    c1=[c1;repmat(yr(i),sum(use1),1)];
    
    gx2=[gx2;dat.gx(use2)];
    gy2=[gy2;dat.gy(use2)];
    dbh2=[dbh2;dat.dbh(use2)];
    agb2=[agb2;dat.agb(use2)];
    ID2=[ID2;dat.treeID(use2)];
    c2=[c2;repmat(yr(i),sum(use2),1)];

    
end
  
% ID0 = unique(ID);
% for i=1:length(ID0)
%     
%     use=ID==ID0(i);
%     plot(c(use),dbh(use),'-k');hold all
% end

B = nan(8,2);
ci1 = nan(8,2);
ci2 = nan(8,2);
for i=2:8
    
    B(i,1)=nansum(agb1(c1==yr(i)))/50*1e-3;
    B(i,2)=nansum(agb2(c2==yr(i)))/50*1e-3;
    
    ci1(i,:) = quantile(bootstrp(1000,'nansum',agb1(c1==yr(i))),[.25 .75])/50*1e-3-B(i,1);
    ci2(i,:) = quantile(bootstrp(1000,'nansum',agb2(c2==yr(i))),[.25 .75])/50*1e-3-B(i,2);
end


clf
errorbar(yr,B(:,1),ci1(:,1),ci1(:,2));hold all
errorbar(yr,B(:,2),ci2(:,1),ci2(:,2))
xlim([1982 2017])
legend('Ceiba pentandra','Dipterix panamensis','location','SouthEast','fontangle','italic')
legend('boxoff')
xlabel('year of cesus')
ylabel('AGB (Mg ha^{-1})')
        
%% load data
clear
dat(1) = load('C:\Users\mdetto\Dropbox (Princeton)\bci_Matlab_Files\bci1.mat');
dat(2) = load('C:\Users\mdetto\Dropbox (Princeton)\bci_Matlab_Files\bci2.mat');
dat(3) = load('C:\Users\mdetto\Dropbox (Princeton)\bci_Matlab_Files\bci3.mat');
dat(4) = load('C:\Users\mdetto\Dropbox (Princeton)\bci_Matlab_Files\bci4.mat');
dat(5) = load('C:\Users\mdetto\Dropbox (Princeton)\bci_Matlab_Files\bci5.mat');
dat(6) = load('C:\Users\mdetto\Dropbox (Princeton)\bci_Matlab_Files\bci6.mat');
dat(7) = load('C:\Users\mdetto\Dropbox (Princeton)\bci_Matlab_Files\bci7.mat');
dat(8) = load('C:\Users\mdetto\Dropbox (Princeton)\bci_Matlab_Files\bci8.mat');
yr = [1982 1985 1990 1995 2000 2005 2010 2015];  
%%  gap=phase dynamics
i0 = 1;
tr = 700;
buf = 1.2;

    use = find(dat(i0).dbh>tr);
    [~,ia,ib] = setxor(dat(i0).treeID(use),dat(i0+1).treeID);
    
    DBH = dat(1).dbh(use(ia));
    SP = dat(1).sp(use(ia));
    clear centers
    centers(:,1) = dat(i0).gx(use(ia));
    centers(:,2) = dat(i0).gy(use(ia));
    radii = buf*sqrt((0.66*(DBH/10).^(1.34))/pi);%crown radius
    
    figure(1)
    clf
    viscircles(centers,radii)
    axis([0 1000 0 500])
    daspect([1 1 1])
    
        
    use = (strcmp(dat(i0+1).sp,'ceibpe') | strcmp(dat(i0+1).sp,'diptpa')) & dat(i0+1).dbh<tr;
    x=dat(i0+1).gx(use);y=dat(i0+1).gy(use);
    hold all
    plot(x,y,'.')
    pause
    
%%     for i=1:length(radii)
b=0;clf
    for i=[27 84 69 83]
        i
        use = (dat(i0).gx-centers(i,1)).^2 + ...
            (dat(i0).gy-centers(i,2)).^2 < (buf*radii(i))^2 & ...
            (strcmp(dat(i0).sp,'ceibpe') | strcmp(dat(i0).sp,'diptpa')) & dat(i0).dbh<DBH(i);
        
        
        if sum(use)>0
            ID1 = dat(i0).treeID(use);
            
            use = (dat(i0).gx-centers(i,1)).^2 + ...
                (dat(i0).gy-centers(i,2)).^2 < (buf*radii(i))^2 & ...
                ~strcmp(dat(i0).sp,'ceibpe') & ~strcmp(dat(i0).sp,'diptpa');
            
            ID2 = dat(i0).treeID(use);
            B1 = nan(length(ID1),8);
            B2 = nan(length(ID2),8);
            for j=i0:8
                
                 [~,ia,ib] = intersect(dat(j).treeID,ID1);
                 B1(ib,j) = dat(j).dbh(ia);
                 
                 [~,ja,jb] = intersect(dat(j).treeID,ID2);
                 B2(jb,j) = dat(j).dbh(ja);
            end
            if i==69;B1(1)=800;end
            
            b=b+1;
            subplot(2,2,b)
            plot(yr,B2/10,'-k');hold all
            z=plot(yr,B1/10,'-r','linewidth',2);hold all
%             h=yline(DBH(i));h.Color='b';
            xlim([1980 2017])
            xlabel('year of cesus')
            ylabel('DBH (cm)')
            set(gca,'ytick',0:50:250)
%             legend(z,dat(j).sp(ia))
           if i==27 || i==84
              legend(z,'\itCeiba pentandra');legend('boxoff');
           else
              legend(z,'\itDipterix oleifera');legend('boxoff');
           end
            pause(.1)
        end
    end
    
    
    
 %% idividual-base analysis
i0 = 1;
tr = 100;
R  = 5;
   
    use = (strcmp(dat(1).sp,'ceibpe') | strcmp(dat(1).sp,'diptpa')) & dat(1).dbh<tr;
    x=dat(1).gx(use);
    y=dat(1).gy(use);
    ID = dat(1).treeID(use);
    
    for i=1:length(ID)
        
        use = ((dat(1).gx-x(i)).^2 + (dat(1).gy-y(i)).^2) < R^2 & ...
            (~strcmp(dat(1).sp,'ceibpe') & ~strcmp(dat(1).sp,'diptpa'));
        ID2 = dat(1).treeID(use);
        
        B1 = nan(8,1);
        B2 = nan(length(ID2),8);
        for j=1:8
            
            [~,ia,ib] = intersect(dat(j).treeID,ID(i));
            if ~isempty(ia)
                B1(j) = dat(j).dbh(ia);
            end
            
            [~,ia,ib] = intersect(dat(j).treeID,ID2);
            if ~isempty(ia)
                B2(ib,j) = dat(j).dbh(ia);
                hhh = B2(ib,j)<B1(j);
                B2(ib(hhh),j)=nan;
            end
        end
        
        
        clf
%         plot(yr,B2,'-k');hold all
%         plot(yr,B1,'-r','linewidth',2);hold all
plot(yr,nansum(0.66*(B2/10).^1.34)/(pi*R^2),'-k');hold all
plot(yr,(0.66*(B1/10).^1.34)/(pi*R^2),'-r')
        xlim([1982 2017])
        pause(.5)
        
    end
    
        
        
        
        

