% This code will read in the text output files from the LAMMPS molecular
% dynamics simulations.
% It then creates a new file of organized data in a textfile in a new
% directory.

% Written by Alex Lohse
% Modified by Noah Paulson 9/18/2014

file = 'dump.deform375_10e7_5aniso_15_edit';
 
fo = fopen( file , 'r');

t_ct = 1;
while t_ct < 1250 && ~feof(fo)
% for t_ct= 1:1000

    fgetl(fo);
    timestep(t_ct) =  fscanf(fo, '%f\n',1);
    fgetl(fo);
    atoms(t_ct) =  fscanf(fo, '%f\n',1);
    fgetl(fo);
    bounds(:,:,t_ct) = fscanf(fo, '%f %f\n',[2 3])';
    fgetl(fo);
    locations(:,:,t_ct) = fscanf(fo, '%d %d %d %f %f %f\n',[6 7999])';

    t_ct = t_ct + 1;   
end

fclose(fo);

save total_info_15.mat timestep atoms bounds locations
