#%% 1
addpath('/Users/nickheaney/Desktop/CoastSnap/Code/GUI');
addpath('/Users/nickheaney/Desktop/CoastSnap-Jupyter-MATLAB')
load('MapShoreline.mat');

#%% 2
P = improfile(xgrid,ygrid,Iplan,transects.x,transects.y);
improfile(xgrid,ygrid,Iplan,transects.x,transects.y)

#%% 3
P(:,:,1)-P(:,:,3);

#%% 4
[pdf_values,pdf_locs] = ksdensity(P(:,:,1)-P(:,:,3));
ksdensity(P(:,:,1)-P(:,:,3))

#%% 5
thresh_weightings = [1/3 2/3];

#%% 6
thresh_otsu = multithresh(P(:,:,1)-P(:,:,3))

#%% 7
I1 = find(pdf_locs<thresh_otsu);
[~,J1] = max(pdf_values(I1));

#%% 8
I2 = find(pdf_locs>thresh_otsu);
[~,J2] = max(pdf_values(I2));

#%% 9
thresh = thresh_weightings(1)*pdf_locs(I1(J1)) + thresh_weightings(2)*pdf_locs(I2(J2))

#%% 10
plot(pdf_locs,pdf_values)
hold on
plot(pdf_locs([I1(J1) I2(J2)]),pdf_values([I1(J1) I2(J2)]),'ro')
YL = ylim;
plot([thresh thresh], YL,'r:','linewidth',2)
plot([thresh_otsu thresh_otsu], YL,'g:','linewidth',2)
xlabel(xlabel_type,'fontsize',10)
ylabel('Counts','fontsize',10)

#%% 11
RminusBdouble = double(Iplan(:,:,1))- double(Iplan(:,:,3));
image(xgrid,ygrid,Iplan)
figure
image(xgrid,ygrid,RminusBdouble)

#%% 12
ROIx = [transects.x(1,:) fliplr(transects.x(2,:))];
ROIy = [transects.y(1,:) fliplr(transects.y(2,:))];
scatter(ROIx, ROIy, 'bx')
axis ij

#%% 13
Imask = ~inpoly([X(:) Y(:)],[ROIx',ROIy']);
RminusBdouble(Imask) = NaN;

#%% 14
c = contours(X,Y,RminusBdouble,[thresh thresh]);
image(xgrid,ygrid,RminusBdouble)
hold on
contour(X,Y,RminusBdouble,[thresh thresh],'r','linewidth',2)
figure
image(xgrid,ygrid,Iplan)
hold on
contour(X,Y,RminusBdouble,[thresh thresh],'r','linewidth',2)

#%% 15
II = find(c(1,:)==thresh);

#%% 16
if II==1 %If only one line
    startI = 2;
    endI = size(c,2);
else
    D = diff(II);
    D = D - 1;
    D = [D (size(c,2) - II(end))];
    [~,J] = max(D); %Select contour that is the longest continuous contour
    if J == 1
        startI = 2;
    else
        startI = 1+J+sum(D(1:J-1));
    end
    endI = startI+D(J)-1;
end
xyz.x = c(1,startI:endI)';
xyz.y = c(2,startI:endI)';
points = [xyz.x xyz.y];

#%% 17
%Now loop through transects to extract shorelines only at the transects
sl.x = NaN(1,length(transects.x));
sl.y = NaN(1,length(transects.y));

#%% 18
angle = atan(diff(transects.y(:,60))/diff(transects.x(:,60)));

#%% 19
%First subtract rotation_center
points_new = points(:,1:2) - repmat([transects.x(1,60) transects.y(1,60)],size(points,1),1);

#%% 20
%Now rotate the points
points_rot = points_new*[cos(angle) -sin(angle); sin(angle) cos(angle)];

#%% 21
plot(points(:,1),points(:,2))
grid on
axis ij
hold on
plot(transects.x(:,60),transects.y(:,60))
A = [transects.x(:,60) transects.y(:,60)];
A_new = A - repmat([transects.x(1,60) transects.y(1,60)],2,1);
A_rot = A_new*[cos(angle) -sin(angle); sin(angle) cos(angle)];
plot(A_rot(:,1),A_rot(:,2))
plot(points_rot(:,1),points_rot(:,2))