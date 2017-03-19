%% ROI from mm to voxel mask
%_________1________
mm1=importdata('/Users/caterinatrainito/Documents/Caterina/statsPMOD_Rat4/mmtable.txt');
mmspm1=mm1(:,[1 3 2]);
mmspm1(:,1)=-mmspm1(:,1);
mmspm1(:,3)=0;
vox1 = floor(mm2vox(mmspm1));

roispm=zeros(size(cestim));
for i=1:size(vox1,1)
    roispm(vox1(i,2),vox1(i,1))=1;
end
roispm=flipud(roispm);
% figure;
% subplot(2,1,1); imagesc(roispm);
% subplot(2,1,2); imagesc(flipud(double(cestim)).*~roispm); colormap('gray');

%_________2________
mm2=importdata('/Users/caterinatrainito/Documents/Caterina/statsPMOD_Rat4/mmtable2.txt');
mmspm2=mm2(:,[1 3 2]);
mmspm2(:,1)=-mmspm2(:,1);
mmspm2(:,3)=0;
vox2 = floor(mm2vox(mmspm2));
roispm2=zeros(size(cestim));
for i=1:size(vox2,1)
    roispm2(vox2(i,2),vox2(i,1))=2;
end
roispm2=flipud(roispm2);
figure;
subplot(2,1,1); imagesc(roispm+roispm2);
subplot(2,1,2); imagesc(flipud(double(cestim)).*~(roispm+roispm2)); colormap('gray');


