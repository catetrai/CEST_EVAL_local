% LOOK AT PET DICOMS
cd('/Users/caterinatrainito/Documents/Caterina/dicoms_PET/1.2.826.0.1.3417726.3.245890.20170116163055356/1.2.826.0.1.3417726.3.132947.20170116163055356');
direc=dir('*.dcm');
[sorted,ix]=sort(cell2mat({direc.datenum}));
dirname={direc.name};
dirname=dirname(ix);
% figure;
for i=1:numel(dirname)
    hdr=dicominfo(dirname{i});
    numim(i,1)=hdr.ImageIndex;
    im(:,:,i)=dicomread(dirname{i});
%     imagesc(double(im))
%     pause;
end


%%%% LOAD PET
pet=analyze75read('/Users/caterinatrainito/Documents/Caterina/dicoms_PET/sMaAm_For_14_06_Glio_dynFmiso_Rat4_10012017-0001-00089-000089.img');
figure;
for i=1:size(pet,3)
    imagesc(double(pet(:,:,i)))
    pause;
end

%%%% TRANSFORMED PET
tform=affine3d(transfmat);
pettrans=double(imwarp(pet,tform));
figure;
for i=1:size(pet,3)
    imagesc(pettrans(:,:,i))
    pause;
end


figure;
for i=1:size(pet,3)
    imagesc(pet(:,:,i)-pettrans(:,:,i))
    pause;
end



%%%% LOAD MR
mr=analyze75read('/Users/caterinatrainito/Documents/Caterina/dicoms_MRanatom/sMaAm_FOR_14_06_Glio_dynFMISO_Rat4_10012017-530001-00001-000001.img');
figure;
for j=1:size(mr,1)
    imagesc(fliplr(rot90(squeeze(double(mr(j,:,:))),-1)))
    pause;
end