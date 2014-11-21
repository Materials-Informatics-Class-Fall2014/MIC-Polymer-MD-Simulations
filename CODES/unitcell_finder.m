%%%
%
%This code is designed to be run after PE_perfect_xtal.m in order to find
%just one unit cell. It can be modified to give the dimensions of a largers
%unit cell. We can then export to Avogadro and PACKMOL when building
%structures.
%
%NOTE: can only be run once after running PE_perfect_xtal.m due to unit conversion.
%
%Alex Lohse
%11/20/2014
%
%%%

%the number of unit cells in each coordinate
uc_x = 3;
uc_y = 1;
uc_z = 1;

%unit cell dimensions
a = 7.42;
b = 4.95;
c = 2.54;

%shifted dimensions due to building the unit cell shifted to begin with
ashift = 0.646;
bshift = 0.561;

allchains = allchains.*10; %convert from nm to Angstroms

cell = (allchains(:,1)<=(uc_x*a+ashift) & allchains(:,2)<=(uc_y*b+bshift)...
    & allchains(:,3)<(uc_z*c));
unitcell=[];
kk = 1;
for ii=1:numel(cell)
    if cell(ii) == 1;
        unitcell(kk,:) = allchains(ii,:);
        kk = kk + 1;
    end
end
clf;
H = plot3(unitcell(:,1), unitcell(:,2), unitcell(:,3));
color = hsv(1);
set(H,'LineStyle','none','Marker','o','MarkerEdgeColor','k','MarkerFaceColor',color(1,:),'MarkerSize',5);
grid on;
axis tight equal;
xlabel x-axis; ylabel y-axis; zlabel z-axis;

%create an xyz file that can be implemented in Avogadro or PACKMOL
file = 'perfect_unitcell.xyz'
fid = fopen(file, 'w');
fprintf(fid,'%d\n',numel(unitcell(:,1)));
fprintf(fid, 'XYZ file input for perfect PE crystal cell\n');

for ii = 1:numel(unitcell(:,1))
    fprintf(fid,'%s %f %f %f\n','C',unitcell(ii,:));
end;
fclose(fid);
