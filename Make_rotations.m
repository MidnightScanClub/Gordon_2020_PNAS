

iterations = 1000;


if isempty(gcp('nocreate'))
    pool = parpool(8);
end

sphere = gifti(['/data/cn/data1/scripts/CIFTI_RELATED/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Conte69.L.sphere.32k_fs_LR.coord.gii']);
sphereLcoords = sphere.vertices;
sphere = gifti(['/data/cn/data1/scripts/CIFTI_RELATED/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Conte69.R.sphere.32k_fs_LR.coord.gii']);
sphereRcoords = sphere.vertices;
allspherecoords = [sphereLcoords ; sphereRcoords];

cifti_template = ft_read_cifti_mod('/data/cn/data1/scripts/CIFTI_RELATED/Resources/Baddata_bigcluster_LR.dtseries.nii');
ncortverts = nnz(cifti_template.brainstructure(cifti_template.brainstructure>0)<3);
nsurfverts = nnz(cifti_template.brainstructure<3);
nsurfvertsL = nsurfverts ./ 2;
nsurfvertsR = nsurfverts ./ 2;

cifti_template.data((ncortverts+1):end) = [];
cifti_template.pos(cifti_template.brainstructure>2,:) = [];
cifti_template.brainstructure(cifti_template.brainstructure>2) = [];
cifti_template.brainstructurelabel(3:end) = [];

origspherecoords = allspherecoords(cifti_template.brainstructure>0 & cifti_template.brainstructure<3,:);
origspherehems = cifti_template.brainstructure(cifti_template.brainstructure>0 & cifti_template.brainstructure<3);
brainstructure_surfacespaceL = cifti_template.brainstructure(1:nsurfvertsL);
brainstructure_surfacespaceR = cifti_template.brainstructure((nsurfvertsL+1):(nsurfvertsL + nsurfvertsR));


remapped = zeros(size(cifti_template.data,1),iterations);

parfor iter = 1:iterations
    
    disp(iter)
    test = zeros(length(brainstructure_surfacespaceL),3);
    
    ok = 0;
    while ok==0
        
        xrot = rand * 2*pi;
        yrot = rand * 2*pi;
        zrot = rand * 2*pi;
        
        if (min([xrot-0, 2*pi - xrot]) + min([yrot-0, 2*pi - yrot]) + min([zrot-0, 2*pi - zrot])) > (pi/8)
            ok = 1;
        end
    end
    
    rotmat_x = [1 0 0;0 cos(xrot) -sin(xrot); 0 sin(xrot) cos(xrot)];
    rotmat_y = [cos(yrot) 0 sin(yrot); 0 1 0; -sin(yrot) 0 cos(yrot)];
    rotmat_z = [cos(zrot) -sin(zrot) 0; sin(zrot) cos(zrot) 0; 0 0 1];
    
    xrotcoords = rotmat_x * origspherecoords';
    xyrotcoords = rotmat_y * xrotcoords;
    xyzrotcoords = rotmat_z * xyrotcoords;
    
    yrot = -yrot;
    zrot = -zrot;
    rotmat_x = [1 0 0;0 cos(xrot) -sin(xrot); 0 sin(xrot) cos(xrot)];
    rotmat_y = [cos(yrot) 0 sin(yrot); 0 1 0; -sin(yrot) 0 cos(yrot)];
    rotmat_z = [cos(zrot) -sin(zrot) 0; sin(zrot) cos(zrot) 0; 0 0 1];
    
    xrotcoords = rotmat_x * origspherecoords';
    xyrotcoords = rotmat_y * xrotcoords;
    xyzrotcoords_mirror = rotmat_z * xyrotcoords;
    
    this_remapped = zeros(size(cifti_template.data,1),1);
    
    for ind = 1:ncortverts
        
        if origspherehems(ind)==1
            test(:,1) = xyzrotcoords(1,ind); test(:,2) = xyzrotcoords(2,ind); test(:,3) = xyzrotcoords(3,ind);
            [~, rot_ind] = min(sum(abs(sphereLcoords - test),2));
            if brainstructure_surfacespaceL(rot_ind) > 0
                rot_cifti_ind = rot_ind - nnz(brainstructure_surfacespaceL(1:rot_ind)<0);
                this_remapped(ind) = rot_cifti_ind;
            end
        else
            test(:,1) = xyzrotcoords_mirror(1,ind); test(:,2) = xyzrotcoords_mirror(2,ind); test(:,3) = xyzrotcoords_mirror(3,ind);
            [~, rot_ind] = min(sum(abs(sphereRcoords - test),2));
            if brainstructure_surfacespaceR(rot_ind) > 0
                rot_cifti_ind = rot_ind - nnz(brainstructure_surfacespaceR(1:rot_ind)<0) + nnz(brainstructure_surfacespaceL>0);
                this_remapped(ind) = rot_cifti_ind;
            end
        end
    end
    remapped(:,iter) = this_remapped;
end

out = cifti_template;
out.data = remapped;
out.dimord = 'pos_time';
ft_write_cifti_mod(['/data/nil-bluearc/GMT/Evan/MSC/reliability_correction/Rotated_inds'],out);


