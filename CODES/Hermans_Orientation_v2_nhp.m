%This is version 2 of the code to calculate Herman's Orientation within a
%voxel. This code integrates Noah Paulson's chain finding algorithm to
%eliminate bad best fit line calculations. The Herman's orientation is
%found by finding the vector between every other monomer in the chain
%(reducing error from bond angle) and comparing that with the average
%vector through all the chains in a given voxel.
%
% Written by Alex Lohse 10/30/2014
% Edited by Noah Paulson 11/02/2014


if ~exist('timestep','var')
    filename = 'tilt6l.mat';
    load(filename)
end

if ~exist( 'assets','dir')
    mkdir('assets');
end

ti = 1; % timestep


data = struct('tme',timestep(ti),...
              'bnds',bounds(:,:,ti),...
              'chids',locations(:,2,ti),...
              'loc',locations(:,4:6,ti));

vloc = ones(size(data.loc(:,1)));

% translate the simulation volume so the lower bottom corner cooincides
% with the 0,0,0 origin
data.loc = data.loc - bsxfun(@times,[vloc,vloc,vloc],data.bnds(:,1,ti)');

val = 1:length(data.loc(:,1)); % IDs for indv. monomers at this timestep

blen = [9,9,9]; % number of Angstroms per bin in the x,y and z directions

% subs: array where each column is a list bin IDs for the chain ID and the
% x, y, z coordinates respectively. For example, the 4th row of subs
% represents the spatial (and chain-wise) binning of monomer 4 where the 
% x position is in bin 2, the y position is in bin 1 and the z position is
% in bin 12
subs = floor(bsxfun(@rdivide,data.loc,blen)) + 1;
subs = [subs, data.chids];

% A: 3d cell array where each cell is identified by an x, y and z location
% as well as a chain ID and contains the IDs of all monomers within that
% binning.
A = accumarray(subs,val,[],@(x) {x});

% next 3 lines: remove cells on the borders of the simulation volume.
A(end,:,:,:) = [];
A(:,end,:,:) = [];
A(:,:,end,:) = [];

color = hsv(20);


%% 3d plot of all chains in different colors

tic

Asize = size(A);
ims = {};

for vox = 1 : (Asize(1)*Asize(2)*Asize(3))

    clf

    voxel_chains = [];
    vectors = [];
    hh = 1;

    % vvec is the spatial location of the current bin in x,y and z
    [vvec(1),vvec(2),vvec(3)] = ind2sub(Asize(1:end-1),vox);

    for chid = 1:20
         
        chains = data.loc(A{vvec(1),vvec(2),vvec(3),chid},:);
        chlen = length(chains(:,1));
        
        figure(1)
        H = plot3(chains(:,1), chains(:,2), chains(:,3));
        set(H,...
            'LineStyle','none',...
            'Marker','o',...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor',color(chid,:),...
            'MarkerSize',5)

        if chid == 1; hold on; end;
               
        if chlen > 2
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

            if isempty(allpair) == 1; break; end; %end loop if no pairs exist
            
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

                markers = ['x','o','s','p','h','^'];
                
                plot3(chains(subchain,1),chains(subchain,2),chains(subchain,3),...
                   'LineStyle','none',...
                   'Marker',markers(xx),...
                   'MarkerSize',10)

                % when there are no more 'ends' left in singles, end the loop,
                % all subchains have been found for the particular chain.
                if isempty(singles) == 1; break; end;

            end
            %This loop generates a variable, voxel_chains, that assigns each of
            %the chains in the voxel a new chain id and then the corresponding
            %location data.
            
            for xx = 1:numel(all_subchain)
                if (numel(all_subchain{xx})) > 2

                    voxel_chains = [voxel_chains;...
                                    hh*ones(numel(all_subchain{xx}),1),...
                                    chains(all_subchain{xx},:)];
                    hh=hh+1;
                end
            end
        end
        
        %empty the subchain array for the next loop to avoid loop counting
        %errors.
        all_subchain = [];
    end

    toc

    if isempty(voxel_chains) == 1; continue; end;
    
    % plot vectors that span every other monomer in a subchain
    for ll = 1:numel(voxel_chains(:,1))-2
        %only do the calculation if the subchain is more than 4 monomers
        if voxel_chains(ll,1)==voxel_chains(ll+2,1)

            coord = [voxel_chains(ll,2:4);voxel_chains(ll+2,2:4)];

            plot3(coord(:,1), coord(:,2), coord(:,3),'k-')

            dirvec = coord(2,:) - coord(1,:);
        end
    end

    % chainmash: precursor step to find vectors in voxel
    chainmash = voxel_chains(3:end,:) - voxel_chains(1:end-2,:);
    % zerpres: vector length of 'voxel_chains' where 1s represent relevant
    % vectors
    zerpres = [chainmash(:,1) == 0];
    % vectors: the set of all direction vectors in the voxel
    vectors = chainmash(zerpres,2:end);
    % negvec: vector length of vectors where 1s represent vectors with a
    % negative z-component
    negvec = [vectors(:,3) < 0];
    % negmat: matrix shape of vectors where rows are negative at positions
    % in negvec with a one
    negmat = ones(length(vectors(:,1)),3) - 2 * [negvec, negvec, negvec];
    % vectors: the set of all direction vectors in the voxel corrected for
    % direction
    vectors = negmat .* vectors; 
    
    
    %calculate the average vector from the matrix of all every other pair
    %vectors
    avgvec = mean(vectors)./norm(mean(vectors));

    %% Plot the mean vector
    middleloc = mean(voxel_chains(:,2:4),1);
    pltavgvec = [middleloc ; middleloc + avgvec];
    plot3(pltavgvec(:,1),pltavgvec(:,2),pltavgvec(:,3),...
        'r-',...
        'linewidth',2)

    
    %% Calculate the Herman's Orientation

    vecnorms = sqrt(sum(vectors.^2,2));

    % the vector angle is calculated by rearanging the formula:
    % cos(theta) = (A dot B) / (||A||*||B||)
    numer = sum(bsxfun(@times,vectors,avgvec),2);
    denom = bsxfun(@times,vecnorms,1);
    angles = abs(acos(numer ./ denom));
    % if the angle is over 90 degrees reverse the vector
    shift = [angles >= 0.5*pi];
    angles = abs(angles - pi*shift);
    
    % Record the set of angles for each voxel
    save_angles{vox} = angles;
    
    % Hermans orientation calc
    Hermans = (3*(cos(angles).^2)-1)/2;

    HermOrientation(vox) = mean(Hermans);

    
    %% Save the plot for each voxel
    hold off
    axis tight equal; grid on;    
    title2 = 'blarg';
    title(title2)   
    
    view(3)
    ims{vox} = fullfile( 'assets', sprintf('chaindir_%i.png',vox) );
    saveas( gcf, ims{vox} )
     
end
