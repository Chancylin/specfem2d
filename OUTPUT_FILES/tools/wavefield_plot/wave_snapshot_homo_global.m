%%load the file which describes the new coordinate of your plotting point
%%because the original output of wavefield is sampled at GLL points what looks 
%%not even distributed. This step is up to you. I just build the points with 
%%constant interval in x and z direction.
coor_new=load('./new_point_global.txt');

%the directory containing wavefield data files
%name2 = '../../global_reconst_full/';
name2 = '../../global_homo_reconst/';
name3 = [name2 'snapshot/'];
%load the coordinate of grids
grid_coor=load([name2 'wavefield_grid_for_dumps_000.txt']);


    
%%configure your color scheme
%%% Enlarge figure to full screen.
% 	set(gcf, 'Units', 'Normalized', 'Outerposition', [0, 0, 1, 1]);
%%I haven't tried too much about different options. You could explore a bit 
%%if you want things look more vivid

%loop for different time step
%for istep = 200:200:6000;
for istep = 1200:600:2400;
    
    figure('Position',[1,1,800,800])
 
     xlim([-10000 10000])
     ylim([-10000 10000])
    
    t=istep*0.005;
    name1=num2str(istep);
    if (istep<1000)
        name1=['00' num2str(istep)];
    elseif (istep<10000)
        name1=['0' num2str(istep)];
    end
    filename=[name2 'wavefield00' name1 '_01_000.txt'];
    field=load(filename);
    %filename_back=[name2_back 'wavefield00' name1 '_01_000.txt'];
    %field_back=load(filename_back);
    %field_scat=field - field_back;
    field_scat=field;
    
%%plot boundaries
%%basically you can ignore the following part because I just draw 
%%those to indicate the geometry of my model.

    

%Interpolate scattered data
    F=TriScatteredInterp(grid_coor(:,1),grid_coor(:,2),field_scat(:,1));
    x_new = F(coor_new(:,1),coor_new(:,2));
%plot the wavefield
    scatter(coor_new(:,1),coor_new(:,2),20,x_new,'filled','o');
    hold on;
    %xlim([-11000 11000])
    %ylim([-11000 11000])
    
    
%%configure your color scheme
    colormap(redblue)
    caxis([-2e-2 2e-2])
    %colorbar;
    %if (istep == 1600)
        colorbar('FontSize',22);%'YTickLabel',{'-2e-2','2e-2'})
    %end
%%% Enlarge figure to full screen.
% 	set(gcf, 'Units', 'Normalized', 'Outerposition', [0, 0, 1, 1]);
%%I haven't tried too much about different options. You could explore a bit 
%%if you want things look more vivid
	caption = sprintf('%.1f s', t);
	%title(caption, 'FontSize', 24);

%%plot boundaries
%%basically you can ignore the following part because I just draw 
%%those to indicate the geometry of my model.
    x = -10000:100:10000;
    y = 10000*ones(1,length(x));
    plot(x,y,'-k','Linewidth',5)
    hold on
    plot(y,x,'-k','Linewidth',5)
    hold on
    x = -10000:100:10000;
    y = -10000*ones(1,length(x));
    plot(x,y,'-k','Linewidth',5)
    hold on
    plot(y,x,'-k','Linewidth',5)
    hold on
    
    x = -5000:100:5000;
    y = 5000*ones(1,length(x));
    plot(x,y,'--k','Linewidth',5)
    hold on
    plot(y,x,'--k','Linewidth',5)
    hold on
    x = -5000:100:5000;
    y = -5000*ones(1,length(x));
    plot(x,y,'--k','Linewidth',5)
    hold on
    plot(y,x,'--k','Linewidth',5)
    hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%drawnow; make the movie
	%frame(k) = getframe(gcf);
    %k=k+1;
    snapshot=[name3 'snapshot' name1];
    set(gca,'XtickLabel',[],'YtickLabel',[]);
    %axis off;
    %set(gca,'visible','off');
    %title(caption, 'FontSize', 24);
    print(snapshot,'-depsc');
    %print(snapshot,'-djpeg');
	% Write this frame out to a new video file.
	%myMovie(frameIndex) = thisFrame;
end
%movie2avi(frame,'wave_movie_true.avi','fps',2);
