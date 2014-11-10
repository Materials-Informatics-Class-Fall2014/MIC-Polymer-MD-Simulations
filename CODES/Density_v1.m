clear;
openfile = '..\Polymer-MD\pe00tilt6l\';
file = 'dump.deform0.1_tilt6l';

fid = fopen(file);
fgetl(fid);
timestep = fscanf(fid, '%d\n');

for ii = 1:3
    fgetl(fid);
end

boundX = fscanf(fid, '%f %f\n', [1 2]);
boundY = fscanf(fid, '%f %f\n', [1 2]);
boundZ = fscanf(fid, '%f %f\n', [1 2]);
bounds = [boundX(1) boundX(2); boundY(1) boundY(2); boundZ(1) boundZ(2)];
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

%stores all xyz data into a variable called locations
locations = fscanf(fid,'%*d %*d %*d %f %f %f\n', [3 7999])';
fclose(fid);

ii = 1;

if ~exist( 'assets','dir')
    mkdir('assets');
end

for ii = 1:number_voxels
    atom_number = (locations(:,1)>voxel_limits{ii}(1,1) & locations(:,1)<voxel_limits{ii}(1,2)...
        & locations(:,2)>voxel_limits{ii}(2,1) & locations(:,2)<voxel_limits{ii}(2,2)...
        & locations(:,3)>voxel_limits{ii}(3,1) & locations(:,3)<voxel_limits{ii}(3,2));
    
    atom_number(atom_number==0) = [];
    
    numin_voxel(ii) = numel(atom_number);
end

perfect_density = 1.004; %data for HDPE, g/cm^3
amorphous_density = 0.853; 
MW_CH2 = 14.027; %molecular weight of ch2, g/mol
Avogadro = 6.022 * 10^23; %avogadros number, atoms/moles
%below volume in Angstrom^3
Volume = (xvoxel(2) - xvoxel(1))*(yvoxel(2)-yvoxel(1))*(zvoxel(2)-zvoxel(1));
Volume = Volume*10^(-24); %in cm^3

for ii = 1:number_voxels
    rho = numin_voxel(ii)*MW_CH2/Avogadro/Volume;
    Density(ii) = rho;
    xtal_degree(ii) = (rho-amorphous_density)/(perfect_density-amorphous_density);
end

NormDens = Density./(max(Density));

%build a polynomial box that we can fill with a single color that will
%represent density
numhsv = 101-100*min(NormDens);
color = hsv(numhsv);
%rgbcolor = hsv2rgb(hsvcolor);
colormap(color);
% set(gca,'CLim',[min(NormDens),1]);

for jj = 1:number_voxels
    xlo = voxel_limits{jj}(1,1);
    xhi = voxel_limits{jj}(1,2);
    ylo = voxel_limits{jj}(2,1);
    yhi = voxel_limits{jj}(2,2);
    zlo = voxel_limits{jj}(3,1);
    zhi = voxel_limits{jj}(3,2);
    
    edges = [xhi - xlo,yhi - ylo,zhi - zlo];
    origin = [xlo,ylo,zlo];
    trans = NormDens(jj);
    index = round(NormDens(jj)*100)-round(min(NormDens)*100)+1;
    col = color(index,:);
    
    for kk = 1:numel(xvoxel)
        if xhi == xvoxel(kk)
%             subplot(numel(xvoxel)-1,1,kk-1);
            plotcube(edges, origin, 1, col);
            axis([xvoxel(kk-1) xvoxel(kk) boundY(1) boundY(2) boundZ(1) boundZ(2)]);
            break;
        end
    end
    hold on;
    
end
colorbar;
axis tight equal;
xlabel x;ylabel y;zlabel z;
view(270,0);
camroll(-90);
hold off;
% for ii = 1:numel(xvoxel)-1
%     subplot(numel(xvoxel),1,ii);
%     view(270,0);
%     camroll(-90);
% end
