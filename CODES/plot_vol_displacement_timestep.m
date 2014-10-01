% This script computes and plots simulation volume and particle
% displacement versus the simulation timestep
%%% Noah Paulson, 9/30/2014

close all

for ii = [1, 2, 3, 4, 5]
    
    filenum = ii + 10;
    file = ['total_info_',int2str(filenum),'.mat'];
    load(file,'timestep','bounds','locations')
    timestep = timestep';
    
    % vol: vector of volumes at each timestep
    vol = squeeze(8*bounds(1,2,:).*bounds(2,2,:).*bounds(3,2,:));  
    
    % these are simply color IDs
    col = [[0,0,0];[0,0,1];[1,0,0];[0,1,1];[0,1,0]];
    
    x_pos = zeros(length(locations(1,1,:)),1);
    y_pos = zeros(length(locations(1,1,:)),1);
    z_pos = zeros(length(locations(1,1,:)),1);
      
    % calculate the positions of a particular monomer though time
    for jj = 1 : length(locations(1,1,:))
        if jj == 1
            x_pos(jj) = 0;
        else
            % atm: the index of the specified monomer
            atm = find(locations(:,1,jj) == 100);
            
            x_pos(jj) = x_pos(jj - 1) + pdist([locations(atm,4,jj);locations(atm,4,jj-1)]);
            y_pos(jj) = y_pos(jj - 1) + pdist([locations(atm,5,jj);locations(atm,5,jj-1)]);
            z_pos(jj) = z_pos(jj - 1) + pdist([locations(atm,6,jj);locations(atm,6,jj-1)]);
        end
    end

    % plot the various figures
    
    figure(1)
    
    h1 = plot(timestep,vol,'Color',col(1,:));
    
    set(gca,'XMinorTick','on','YMinorTick','on')
    title1 = ['Simulation volume versus timestep, aniso: ', int2str(filenum)];
    title(title1)
    xlabel 'timestep (fs)'
    ylabel 'simulation volume'
    axis([min(timestep) max(timestep) 2E5 5E5])
    grid on
    
    file = ['volume_timestep_aniso', int2str(filenum), '.png'];
    saveas(h1, file)

    figure(2)
    
    plot(timestep,x_pos,'Color',col(1,:))
    hold on
    plot(timestep,y_pos,'Color',col(2,:))
    h2 = plot(timestep,z_pos,'Color',col(3,:));
    hold off
    
    legend('x-displacement','y-displacement','z-displacement')
    title2 = ['Displacements for single monomer versus timestep, aniso: ', int2str(filenum)];
    title(title2)
    xlabel 'timestep (fs)'
    ylabel 'displacements'
    axis([min(timestep) max(timestep) 0 10000])
    grid on

    file = ['displacement_timestep_aniso', int2str(filenum), '.png'];
    saveas(h2, file)
    
    
    figure(3)
    
    plot(x_pos,vol,'Color',col(1,:))
    hold on
    plot(y_pos,vol,'Color',col(2,:))
    h3 = plot(z_pos,vol,'Color',col(3,:));
    hold off
    
    legend('x-displacement','y-displacement','z-displacement')
    title3 = ['Displacements versus simulation volume, aniso: ', int2str(filenum)];
    title(title3)
    xlabel 'displacements'
    ylabel 'simulation volume'
    axis([0 10000 2E5 5E5])
    grid on
    
    file = ['displacement_volume_aniso', int2str(filenum), '.png'];
    saveas(h3, file)
    
end
