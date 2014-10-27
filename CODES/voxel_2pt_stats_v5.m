if exist('timestep','var') == 0
    load tilt6l.mat
end

close all
color = hsv(20);
ti = 1;

% lines which plot a box centered at 0,0,0 with side lengths of 1 unit
box =   [ .5, .5,-.5;...
         -.5, .5,-.5;...
         -.5,-.5,-.5;...
          .5,-.5,-.5;...
          .5, .5,-.5;...
          .5, .5, .5;...
         -.5, .5, .5;...
         -.5,-.5, .5;...
          .5,-.5, .5;...
          .5, .5, .5;...
         -.5, .5, .5;...
         -.5, .5,-.5;...
         -.5, .5, .5;...
         -.5,-.5, .5;...
         -.5,-.5,-.5;... 
         -.5,-.5, .5;...
          .5,-.5, .5;...
          .5,-.5,-.5;...
          .5,-.5, .5]; 

blen = 9;
xdim = floor((2*bounds(1,2,ti))/blen);
ydim = xdim;
zdim = floor((2*bounds(3,2,ti))/blen);

loc = locations(:,4:6,ti);

indc = zeros((xdim+1)*(ydim+1)*(zdim+1),3);

k = 1;

vloc = ones(size(loc(:,1)));
translation = [bounds(1,1,ti)*vloc,bounds(2,1,ti)*vloc,bounds(3,1,ti)*vloc];
loc = loc - translation;
cmax =(xdim+1)*(ydim+1)*(zdim+1);


for chid = 1:20

    chain = (locations(:,2,ti) == chid);
    xvox = chain .* loc(:,1);
    yvox = chain .* loc(:,2);
    zvox = chain .* loc(:,3);
    xvox(xvox==0) = [];
    yvox(yvox==0) = [];
    zvox(zvox==0) = [];
    
    chains = [xvox,yvox,zvox];
    
    % 3d plot of all chains
    figure(1)
    
    H = plot3(chains(:,3), chains(:,1), chains(:,2));
    set(H,'LineStyle','none','Marker','o','MarkerEdgeColor','k',...
        'MarkerFaceColor',color(chid,:),'MarkerSize',5)
    axis tight equal; grid on;    
    title2 = ['All Chains, Timestep = ', num2str(timestep(ti))];
    title(title2)

    if chid == 1; hold on; end;
    
end

tic
c = 1;
for xx = 0:xdim-1
    
    lx = xx * blen;
    xpos = (loc(:,1) > lx & loc(:,1) < (lx + blen));
    
    for yy = 0:ydim-1
        
        ly = yy * blen;
        ypos = (loc(:,2) > ly & loc(:,2) < (ly + blen));
        
        for zz = 0:zdim-1
        
            lz = zz * blen;
            zpos = (loc(:,3) > lz & loc(:,3) < (lz + blen));

            % vector of 1 and 0s where the one designates the preceding
            % conditions have been met. This will be multiplied by the
            % coordinates of the monomers to only get a monomer from  a
            % certain chain in a certain spatial bin.
            mult0 = xpos .* ypos .* zpos;

            xvox = mult0 .* loc(:,1);
            yvox = mult0 .* loc(:,2);
            zvox = mult0 .* loc(:,3);
            xvox(xvox==0) = [];
            yvox(yvox==0) = [];
            zvox(zvox==0) = []; 
            
            vox = [xvox,yvox,zvox];
            all_vox{c} = vox;
            
            indc(c,:) = [xx,yy,zz];
            c = c + 1;
           
            % Plot boxes around the 8 corners of the voxelized simulation
            % volume
            if (xx == 0 || xx == (xdim-1)) && ...
               (yy == 0 || yy == (ydim-1)) && ...
               (zz == 0 || zz == (zdim-1))
                
                xs = blen*box(:,1) + lx + 0.5*blen;
                ys = blen*box(:,2) + ly + 0.5*blen;
                zs = blen*box(:,3) + lz + 0.5*blen;
                
                plot3(zs,xs,ys,...
                    'k-','LineWidth',2.0)
                axis tight equal; grid on;     
                hold on
            end   
        end
    end
end
hold off

toc

vid = 80;
scale = 4;

vox0 = all_vox{vid};
voxv = ones(size(vox0(:,1)));
translation = [indc(vid,1)*blen*voxv,...
               indc(vid,2)*blen*voxv,...
               indc(vid,3)*blen*voxv];
vox0b = vox0 - translation;
vox0c = scale*vox0b + 1;
vox0d = round(vox0c);
el = blen*scale + 1;

ms = zeros(el,el,el);
for aa = 1 : length(voxv)
    ms(vox0d(aa,1),vox0d(aa,2),vox0d(aa,3)) = 1;
    
    % Plot the current polymer in the voxel
    figure(2)
    
    
    subplot(1,2,1)
    plot3(vox0c(aa,1),vox0c(aa,2),vox0c(aa,3),...
        'LineStyle','none','Marker','o','MarkerEdgeColor','k',...
        'MarkerFaceColor',color(10,:),'MarkerSize',5);
    title('Placement of Monomer Position within Subvoxels')
    hold on
    
    subplot(1,2,2)
    plot3(vox0c(aa,1),vox0c(aa,2),vox0c(aa,3),...
        'LineStyle','none','Marker','o','MarkerEdgeColor','k',...
        'MarkerFaceColor',color(10,:),'MarkerSize',5);
    hold on
    
    % plot the box which represents the closest bin to each particle
    xs = box(:,1) + vox0d(aa,1);
    ys = box(:,2) + vox0d(aa,2);
    zs = box(:,3) + vox0d(aa,3);
    
    subplot(1,2,1)
    plot3(xs,ys,zs,...
        'b-','LineWidth',0.5)
    axis equal; grid on;
    axis([0.5 (blen*scale + 0.5)...
          0.5 (blen*scale + 0.5)...
          0.5 (blen*scale + 0.5)])
    hold on
    
    subplot(1,2,2)
    plot3(xs,ys,zs,...
        'b-','LineWidth',0.5)
    axis equal; grid on;
    axis([0.5 (blen*scale + 0.5)...
          0.5 (blen*scale + 0.5)...
          0.5 (blen*scale + 0.5)])
    hold on
end
hold off

% Calculate Autocorrelation of particles
[T, xx] = SpatialStatsFFT(ms,[],'display',false,'shift',true);

% Plot point where 2pt stats are non-zero
Timg = find(T > 0.0000001);
sc = zeros(length(Timg),3);
[sc(:,1), sc(:,2), sc(:,3)] = ind2sub(size(T), Timg);

figure(3)

subplot(1,2,1)
plot3(sc(:,1),sc(:,2),sc(:,3),...
    'LineStyle','none','Marker','o','MarkerEdgeColor','k',...
    'MarkerFaceColor',color(17,:),'MarkerSize',5);
axis tight equal;
title('2-point statistics for individual voxel')

subplot(1,2,2)
plot3(sc(:,1),sc(:,2),sc(:,3),...
    'LineStyle','none','Marker','o','MarkerEdgeColor','k',...
    'MarkerFaceColor',color(17,:),'MarkerSize',5);
axis tight equal;


% Plot slice of 2-pt stats as image
figure(4)
image(T(:,:,ceil(0.5*el)),'CDataMapping','scaled')
colormap('jet')
axis tight equal;
colorbar
shading flat
caxis([0 1E-5])
