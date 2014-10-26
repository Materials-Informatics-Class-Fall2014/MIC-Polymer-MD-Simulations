%This code develops perfect crystalline structure based on polyethylene
%information gathered from the Handbook of Polyethylene.

%For reference: PE forms a BC orthrhombic crystal structure

bond_angle = deg2rad(112);
bond_length = 0.153; %in nanometers

%initial dimensions of cell in simulation (nm)
xdim = 3;
ydim = 3;
zdim = 5;

%create a vector with the initial chain
%this vector can then be repeated at the unit cell corners to create a
%complete crystal over the volume
k = 1;
for ii = 0:.127:zdim
    if rem(ii,.254)==0
        initchain(k,1) = 0; %x location
        initchain(k,2) = 0; %y location
        initchain(k,3) = ii; %z location
        bc_chain(k,1) = .4356;
        bc_chain(k,2) = .2475;
        bc_chain(k,3) = ii;
    else
        initchain(k,1) = 0.0646;
        initchain(k,2) = 0.0561;
        initchain(k,3) = ii;
        bc_chain(k,1) = .371;
        bc_chain(k,2) = .3037;
        bc_chain(k,3) = ii;
    end
    k=k+1;
end

%for jj = 0.127:0.254:zdim
%    initchain(k,1) = 0.0646;
%    initchain(k,2) = 0.0561;
%    initchain(k,3) = jj;
%    bc_chain(k,1) = .371;
%    bc_chain(k,2) = .3037;
%    bc_chain(k,3) = ii;
%    k=k+1;
%end

%we now have one complete chain so we can repeat in a and b accordingly
%a rotated chain will need to be created for the (body-center)
chains = vertcat(initchain,bc_chain);
ctr = 1;
for ii = 0.742:0.742:xdim
    newchainx = initchain(:,1)+ctr*0.742;
    newbc_chainx = bc_chain(:,1)+ctr*.742;
    chains = vertcat(chains,[newchainx,initchain(:,2),initchain(:,3)],[newbc_chainx,bc_chain(:,2),bc_chain(:,3)]);
    ctr = ctr + 1;
end;
ctr = 1;
allchains = chains;
for ii = .495:.495:ydim
    newchainy = chains(:,2)+ctr*0.495;
    allchains = vertcat(allchains,[chains(:,1),newchainy,chains(:,3)]);
    ctr = ctr+1;
end
j = 1;
for ii = 1:numel(allchains(:,1))
    atoms(ii) = ii;
    id(ii) = j;
    if rem(ii,numel(initchain(:,1)))==0
        j = j+1; %will create unique chain id's
    end
    mol(ii) = 1;
end
num = numel(initchain(:,1));
col = numel(allchains(:,1))/num;
color = hsv(col);
H = plot3(allchains(1:num,1),allchains(1:num,2),allchains(1:num,3));
set(H,'LineStyle','none','Marker','o','MarkerEdgeColor','k','MarkerFaceColor',color(1,:),'MarkerSize',5);
grid on;
hold on;
for m = 2:col
    H = plot3(allchains((m-1)*num+1:m*num,1),allchains((m-1)*num+1:m*num,2),allchains((m-1)*num+1:m*num,3));
    set(H,'LineStyle','none','Marker','o','MarkerEdgeColor','k','MarkerFaceColor',color(m,:),'MarkerSize',5);
end
axis tight equal;

hold off;
%fid = fopen('perfect-xtal.txt','w');
%fprintf(fid,'ITEM: TIMESTEP\n6400000\n');
%fprintf(fid,'ITEM: NUMBER OF ATOMS\n7999\n');
%fprintf(fid,'ITEM: BOX BOUNDS\n');
%box_bounds = [0.0,4.4;0.0,4.4;0.0,11.2];
%fprintf(fid,'%.6f %.6f\n', box_bounds(:,1), box_bounds(:,2));
%fprintf(fid,'ITEM: ATOMS id mol type x y z \n');
%for ii = 1:numel(atoms)
%    fprintf(fid,'%d %d %d %.6f %.6f %.6f\n',atoms(ii),id(ii),mol(ii),allchains(ii,1),allchains(ii,2),allchains(ii,3));
%end
%fclose(fid);



