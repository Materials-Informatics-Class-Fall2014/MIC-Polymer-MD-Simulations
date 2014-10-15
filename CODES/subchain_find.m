%%%
% This code takes a single time-step from a PE-MD simulatation, generates a
% single spatial bin, plots the chains in the bin, and finally identifies
% which sections of chain in the bin are 'whole' within that bin (ie
% sections not cut by the boundary of the spatial bin)
%
% Noah Paulson 10/14/2014
%%%


close('all')

if exist('timestep') == 0
    load total_info_15.mat
end
    
color = hsv(20);

ii = 1;

for jj = 1:20

    % starting origin for bin of interest
    lx = 0; ly = 0; lz = 0;

    % square dimension of bin
    range = 10;

    chain = (locations(:,2,ii) == jj);
    xpos = (locations(:,4,ii) > lx & locations(:,4,ii) < (lx + range));
    ypos = (locations(:,5,ii) > ly & locations(:,5,ii) < (ly + range));
    zpos = (locations(:,6,ii) > lz & locations(:,6,ii) < (lz + range));

    % vector of 1 and 0s where the one designates the preceding
    % conditions have been met. This will be multiplied by the
    % coordinates of the monomers to only get a monomer from  a
    % certain chain in a certain spatial bin.
    mult0 = chain .* xpos .* ypos .* zpos;

    x_chain = mult0 .* locations(:,4,ii);
    y_chain = mult0 .* locations(:,5,ii);
    z_chain = mult0 .* locations(:,6,ii);
    x_chain(x_chain==0) = [];
    y_chain(y_chain==0) = [];
    z_chain(z_chain==0) = [];

    % plot the actual chains in the bin
    figure(1)
    H = plot3(x_chain, y_chain, z_chain);
    set(H,'LineStyle','none','Marker','o','MarkerEdgeColor','k','MarkerFaceColor',color(jj,:),'MarkerSize',5)
    axis tight equal;
    axis ([lx (lx + range) ly (ly + range) lz (lz + range)])
    grid on;    
    title(horzcat('Timestep = ', num2str(timestep(ii))));
    xlabel x-axis; ylabel y-axis; zlabel z-axis;

    if jj == 1
        hold on;
    end

    chlen = length(x_chain);
    
    chains = [x_chain,y_chain,z_chain];
    
    if chlen > 4

        % find the distances between the particles in a single chain
        dists = pdist(chains); dists = squareform(dists);

        %% generate a list of pairs of monomers within the chain
        allpair = [];
        for kk = 1 : chlen - 1

            % find the pairings where the distance between particles is
            % less than 2 units
            pair_ids = find(dists(kk,(kk+1):end) < 1.65);
            pair_ids = pair_ids + kk * ones(size(pair_ids));

            % add each of these pairings to an array for storage
            for ll = 1 : length(pair_ids)
                allpair = [allpair ; kk, pair_ids(ll)]; 
            end               
        end

        %% find the id's of monomers which are at the end of a chain
        singles = [];
        for kk = 1 : max(allpair(:))

            % the first indices of pairs containing a monomer of id == kk
            currindx = find(allpair(:,1) == kk | allpair(:,2) == kk);
            % the second indices associated with a monomer of id == kk
            position = find(allpair(currindx,:) == kk);                

            % append kk to singles if it is a subchain ending
            if length(currindx) == 1
                singles = [singles; kk];
            end

        end
        
        %% find the subchains
        for xx = 1:100

            prev = singles(1); % pick the first chain start available
            subchain = []; % initialize a subchain
            allpair_red = allpair;

            for kk = 1:100

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
            figure(1)
            markers = ['x','o','s'];
            plot3(x_chain(subchain),y_chain(subchain),z_chain(subchain),...
                'LineStyle','none','Marker',markers(xx),'MarkerSize',10)
            
            pause(4)
            
            % when there are no more 'ends' left in singles, end the loop,
            % all subchains have been found for the particular chain.
            if isempty(singles) == 1; break; end;
            
        end
    end
end

hold off;
