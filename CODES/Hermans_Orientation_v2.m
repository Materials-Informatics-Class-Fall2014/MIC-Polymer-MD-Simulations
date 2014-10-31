%This is version 2 of the code to calculate Herman's Orientation within a
%voxel. This code integrates Noah Paulson's chain finding algorithm to
%eliminate bad best fit line calculations. The Herman's orientation is
%found by finding the vector between every other monomer in the chain
%(reducing error from bond angle) and comparing that with the average
%vector through all the chains in a given voxel.
%
%Alex Lohse 10/30/2014
%

clear;
openfile = '..\Polymer-MD\pe00tilt6l\';
file = 'dump.deform0.1_tilt6l';

fid = fopen(file);
fgetl(fid);
time = fscanf(fid, '%d\n');

for ii = 1:3
    fgetl(fid);
end

boundX = fscanf(fid, '%f %f\n', [1 2]);
boundY = fscanf(fid, '%f %f\n', [1 2]);
boundZ = fscanf(fid, '%f %f\n', [1 2]);
xvoxel = linspace(boundX(1),boundX(2), 4);
yvoxel = linspace(boundY(1),boundY(2), 4);
zvoxel = linspace(boundZ(1),boundZ(2), 20);

%calculates the number of voxels
number_voxels = ((numel(xvoxel)-1)*(numel(yvoxel)-1)*(numel(zvoxel)-1));

voxel_limits = cell(number_voxels,1);
jj = 1;
for j = 1:(numel(xvoxel)-1)
    for k = 1:(numel(yvoxel)-1)
        for l = 1:(numel(zvoxel)-1)
            voxel_limits{jj} = [xvoxel(j),xvoxel(j+1);...
                yvoxel(k),yvoxel(k+1);...
                zvoxel(l),zvoxel(l+1)];
            jj = jj + 1;
        end
    end
end

fgetl(fid);

%stores all xyz data into a variable called data
locations = fscanf(fid,'%d %d %d %f %f %f\n', [6 7999])';
fclose(fid);
color = hsv(20);
ii = 1;
for jj = 1:number_voxels
    %sets the limits of the voxel currently working in
    xlo = voxel_limits{jj}(1,1);
    xhi = voxel_limits{jj}(1,2);
    ylo = voxel_limits{jj}(2,1);
    yhi = voxel_limits{jj}(2,2);
    zlo = voxel_limits{jj}(3,1);
    zhi = voxel_limits{jj}(3,2);
    hh=1;
    voxel_chains = [];
    startpoint = [];
    vectors = [];
    for kk = 1:20
        chain = (locations(:,2) == kk);
        xpos = (locations(:,4) > xlo & locations(:,4) < xhi);
        ypos = (locations(:,5) > ylo & locations(:,5) < yhi);
        zpos = (locations(:,6) > zlo & locations(:,6) < zhi);
        
        % vector of 1 and 0s where the one designates the preceding
        % conditions have been met. This will be multiplied by the
        % coordinates of the monomers to only get a monomer from  a
        % certain chain in a certain spatial bin.
        mult0 = chain .* xpos .* ypos .* zpos;
        
        x_chain = mult0 .* locations(:,4);
        y_chain = mult0 .* locations(:,5);
        z_chain = mult0 .* locations(:,6);
        x_chain(x_chain==0) = [];
        y_chain(y_chain==0) = [];
        z_chain(z_chain==0) = [];
        
        % plot the actual chains in the bin
        %figure(1)
        %H = plot3(x_chain, y_chain, z_chain);
        %set(H,'LineStyle','none','Marker','o','MarkerEdgeColor','k','MarkerFaceColor',color(jj,:),'MarkerSize',5)
        %axis tight equal;
        %axis ([xlo xhi ylo yhi zlo zhi])
        %grid on;
        %title(horzcat('Timestep = ', num2str(time)));
        %xlabel x-axis; ylabel y-axis; zlabel z-axis;
        
        %if jj == 1
        %    hold on;
        %end
        
        chlen = length(x_chain);
    
        chains = [x_chain,y_chain,z_chain];
        
        if chlen > 4
            % find the distances between the particles in a single chain
            dists = pdist(chains); dists = squareform(dists);
            
            %% generate a list of pairs of monomers within the chain
            allpair = [];
            for ll = 1 : chlen - 1
                
                % find the pairings where the distance between particles is
                % less than 2 units
                pair_ids = find(dists(ll,(ll+1):end) < 1.65);
                pair_ids = pair_ids + ll * ones(size(pair_ids));
                
                % add each of these pairings to an array for storage
                for mm = 1 : length(pair_ids)
                    allpair = [allpair ; ll, pair_ids(mm)];
                end
            end
            
            %% find the id's of monomers which are at the end of a chain
            singles = [];
            for nn = 1 : max(allpair(:))
                
                % the first indices of pairs containing a monomer of id ==
                % nn
                currindx = find(allpair(:,1) == nn | allpair(:,2) == nn);
                % the second indices associated with a monomer of id == nn
                position = find(allpair(currindx,:) == nn);
                
                % append kk to singles if it is a subchain ending
                if length(currindx) == 1
                    singles = [singles; nn];
                end
                
            end
            
            %% find the subchains
            for xx = 1:100
                
                prev = singles(1); % pick the first chain start available
                subchain = []; % initialize a subchain
                allpair_red = allpair;
                
                for ll = 1:100
                    
                    % append the previous monomer id to the subchain (in the
                    % first iteration this is the chain-start)
                    subchain = [subchain; prev];
                    
                    % find the first index of the pair associated with the
                    % previous monomer id
                    currindx = find(allpair_red(:,1) == prev | allpair_red(:,2) == prev);
                    % find the second index of the pair of the previous monomer
                    position = find(allpair_red(currindx,:) ~= prev);
                    
                    % if the monomer is not present in the current list of
                    % monomers, end the loop. This means a subchain has been
                    % found
                    if isempty(currindx) == 1; break; end;
                    
                    % the monomer linked to the previous monomer
                    curr = allpair_red(currindx,position);
                    
                    % remove the pairs associated with the previous monomer in
                    % the chain
                    other_indices = find(allpair_red(:,1) ~= prev & allpair_red(:,2) ~= prev);
                    allpair_red = allpair_red(other_indices,:);
                    
                    % assign the id of the current monomer in the chain to the
                    % variable for the previous monomer
                    prev = curr;
                    
                end
                
                % add the subchain to a list of subchain
                all_subchain{xx} = subchain;
                
                % remove the start and end of the previous subchain from
                % singles
                old_indices = find(singles ~= subchain(1) & singles ~= subchain(end));
                singles = singles(old_indices);
                
                % plot markers on top of the original chain to show the unique
                % subchains
                %figure(1)
                %markers = ['x','o','s'];
                %plot3(x_chain(subchain),y_chain(subchain),z_chain(subchain),...
                %    'LineStyle','none','Marker',markers(xx),'MarkerSize',10)
                
                pause(4)
                
                % when there are no more 'ends' left in singles, end the loop,
                % all subchains have been found for the particular chain.
                if isempty(singles) == 1; break; end;
                
            end
            %This loop generates a variable, voxel_chains, that assigns each of
            %the chains in the voxel a new chain id and then the corresponding
            %location data.
            
            for oo = 1:numel(all_subchain)
                if (numel(all_subchain{oo}))>4
                    for pp=1:numel(all_subchain{oo})
                            locdata = chains(all_subchain{oo}(pp),:);
                            voxel_chains = vertcat(voxel_chains,[hh,locdata(1),locdata(2),locdata(3)]);       
                    end
                    hh=hh+1;
                end
            end
        end
        
        %empty the subchain array for the next loop to avoid loop counting
        %errors.
        all_subchain = [];
    end
    
    %Iterate through the voxel chains and calculate the vector between
    %every other pair
    for ll = 1:numel(voxel_chains(:,1))-2
        %only do the calculation if the subchain is more than 4 monomers
        if voxel_chains(ll,1)==voxel_chains(ll+2,1)
            x = [voxel_chains(ll,2);voxel_chains(ll+2,2)];
            y = [voxel_chains(ll,3);voxel_chains(ll+2,3)];
            z = [voxel_chains(ll,4);voxel_chains(ll+2,4)];
            
            [n,mx] = size(x); [ny,my] = size(y); [nz,mz] = size(z);
            %m is the startpoint
            m = [mean(x), mean(y), mean(z)];
            %w is the orthogonal distance between the average
            w = [x-m(1), y-m(2), z-m(3)];
            a = (1/n)*w'*w;
            [u,d,v] = svd(a); 
            %q is the direction vector
            q = u(:,1); %extracts the largest eigenvector as direction
            startpoint = vertcat(startpoint,m);
            vectors = vertcat(vectors,q');
        end
    
    end
    %calculate the average vector from the matrix of all every other pair
    %vectors
    avgstart = mean(startpoint);
    avgvector = mean(vectors);
    
    for ll = 1:numel(vectors(:,1))
        angle = atan2(norm(cross(vectors(ll,:),avgvector)),dot(vectors(ll,:),avgvector));
        %if angle > (pi/2)
         %   angle = pi - angle;
        %end
        AllAngles(ll) = angle;
    end
    
    Hermans = (3*(cos(AllAngles(:)).^2)-1)/2;
    
    HermOrientation(jj) = mean(Hermans)
    
end
    
%hold off;
