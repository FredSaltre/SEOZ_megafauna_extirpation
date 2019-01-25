program Megafauna_analyses

%%=========================================================================
%% STEP #1: MEGAFAUNA INPUTS
%%=========================================================================  
clear all, close all, clc
new=load('AllMegafauna_Finalcutoff_v2.txt');old=load('AllMegafauna_Finalcutoff_v1.txt');
seaclr=brewermap(15,'RdBu');seaclr1=seaclr(9,:);clr=brewermap(12,'Paired');col=[clr(7,:);clr(4,:)];
coast = load('coast');seaclr=brewermap(15,'RdBu');seaclr1=seaclr(9,:);
figure('Color','w'),worldmap('Australia');gridm('on');setm(gca,'ParallelLabel','off','MeridianLabel','off');
                    hold on,geoshow(flipud(coast.lat),flipud(coast.long),'DisplayType','polygon','FaceColor',seaclr1),
                    hold on,plotm(new(:,2),new(:,1),'ok','Markersize',5,'Markerfacecolor',clr(5,:),'Markeredgecolor',clr(6,:));       
                    hold on,plotm(old(:,2),old(:,1),'pk','Markersize',5,'Markerfacecolor',clr(1,:),'Markeredgecolor',clr(2,:));       
                    tr1=title('Comparison timing Human & Megafauna');set(tr1,'Fontsize',13);legend('','','','','test','test2','location','southwestoutside')
                    
figure('Color','w'),worldmap('Australia');gridm('on');setm(gca,'ParallelLabel','off','MeridianLabel','off');
                    hold on,geoshow(flipud(coast.lat),flipud(coast.long),'DisplayType','polygon','FaceColor',seaclr1),
                    hold on,scatterm(new(:,2),new(:,1),50,new(:,3),'filled');   
                    
%%=========================================================================
%% STEP #2: MEGAFAUNA TIMING OUTPUTS
%%=========================================================================  
clear all, close all, clc,
mat=load('Outputs_AllMegafauna_brut_v3.txt');mask=load('SAHULmask(1x1).txt');nm=length(mat(:,1));xg=-179.5:1:179.5;yg=-89.5:89.5;clat=flipud(yg');clon=xg;
[XI,YI]=meshgrid(xg,yg);YI=flipud(YI);mres=griddata(mat(:,1),mat(:,2),mat(:,3),XI,YI,'v4');sres=griddata(mat(:,1),mat(:,2),mat(:,4),XI,YI,'v4');%interpolation des resultats moyens (mres) + barres d'erreurs (sres)
mres=mres.*mask;[row,col]=find(mres<0);for i=1:length(row),mres(row(i),col(i))=NaN;end;clear row col;%masque les oceans tenant compte du niveau marin de l'epoque et supprime les valeurs negatives
sres=sres.*mask;[row,col]=find(sres<0);for i=1:length(row),sres(row(i),col(i))=NaN;end;clear row col;
TF=isnan(mask);[row,col]=find(mask==1);out=nan(length(row),7);
for i=1:length(row),out(i,:)=[YI(row(i),col(i)),XI(row(i),col(i)),row(i),col(i),mres(row(i),col(i))-sres(row(i),col(i)),mres(row(i),col(i)),mres(row(i),col(i))+sres(row(i),col(i))];end;%[,1]=lat,[,2]=lon,[,3]=row,[,4]=col,[,5]=mean-sd,[,6]=mean,[,7]=mean+sd
save -ascii Megafauna_Fullestimate_v3.txt out;

%% MEGAFAUNA TIMING  ********************************************************** 
clear all, close all, clc,
mat=load('Outputs_AllMegafauna_brut_v3.txt');mask=load('SAHULmask(1x1).txt');nm=length(mat(:,1));xg=-179.5:1:179.5;yg=-89.5:89.5;clat=flipud(yg');clon=xg;
[XI,YI]=meshgrid(xg,yg);YI=flipud(YI);mres=griddata(mat(:,1),mat(:,2),mat(:,3),XI,YI,'v4');sres=griddata(mat(:,1),mat(:,2),mat(:,4),XI,YI,'v4');%interpolation des resultats moyens (mres) + barres d'erreurs (sres)
[FX,FY] = gradient(mres); %identification of spatial gradients for display
mres=mres.*mask;[row,col]=find(mres<0);for i=1:length(row),mres(row(i),col(i))=NaN;end;clear row col;%masque les oceans tenant compte du niveau marin de l'epoque et supprime les valeurs negatives
 
cutlat=find(clat==-26.5);cutlon=find(clon==134.5);mres(1:cutlat,:)=NaN;mres(:,1:cutlon)=NaN;mres(:,340:end)=NaN;save -ascii AllMegafaunaEstimate_v3.txt mres;
dat=load('AllMegafauna_Finalcutoff_v2.txt');clr=brewermap(12,'Paired');clr2=brewermap(13,'Greys');%loading dataset
%glat=load('Latgrad_MegaExtinction(raster1x1).txt');glon=load('Longrad_MegaExtinction(raster1x1).txt');% affichage des gradients latitudinaux and longitudinaux
TF=isnan(mres);[r2,c2]=find(TF==0);nr2=length(r2);glatlon=nan(nr2,4);for i=1:nr2;glatlon(i,:)=[clat(r2(i)),clon(c2(i)),FY(r2(i),c2(i)),FX(r2(i),c2(i))];end;
per=round((length(glatlon(:,1)).*60)/100);pk=randi(nr2,per,1);glatlon2=glatlon(pk,:);
coast = load('coast');seaclr=brewermap(15,'RdBu');seaclr1=seaclr(9,:);col=brewermap(13,'YlOrRd');col=flipud(col);col(1,:)=[];col(end-1:end,:)=[];
figure('Color','w'),worldmap('Australia');gridm('on');setm(gca,'ParallelLabel','off','MeridianLabel','off');
                    hold on,geoshow(flipud(coast.lat),flipud(coast.long),'DisplayType','polygon','FaceColor',seaclr1),
                    hold on,surfm(clat,clon,mres);hold on,plotm(coast.lat,coast.long,'-k');colormap(col);caxis([30000 50000]);
                    hold on,ar=quiverm(glatlon2(:,1),glatlon2(:,2),glatlon2(:,3),-glatlon2(:,4),1.5),set(ar,'color',clr2(7,:),'linewidth',0.5)
                    hold on,plotm(dat(:,2),dat(:,1),'ok','Markersize',3,'Markerfacecolor',clr(1,:),'Markeredgecolor',clr(2,:)); 
                    tr1=title('timing of regional Megafauna extirpation');set(tr1,'Fontsize',13);legend('','','','','','','test2','location','southwestoutside')
                    
%%=========================================================================
%% STEP #3: MEGAFAUNA BEARING OUTPUTS
%%=========================================================================  
%% Calculation of spatial gradients and bearing
clear all, close all, clc,
szlon=deg2km(1.0028,'earth');wrdmulti=load('World_areamultiplier.txt');
mat=load('Outputs_AllMegafauna_brut_v3.txt');mask=load('SAHULmask(1x1).txt');xg=-179.5:1:179.5;yg=-89.5:89.5;
[XI,YI]=meshgrid(xg,yg);YI=flipud(YI);mres=griddata(mat(:,1),mat(:,2),mat(:,3),XI,YI,'v4');%interpolation des resultats moyens (mres) 
mres=mres.*mask;[row,col]=find(mres<0);for i=1:length(row),mres(row(i),col(i))=NaN;end;clear row col;%masque les oceans tenant compte du niveau marin de l'epoque et supprime les valeurs negatives
[n,c]=size(mres);out=zeros(n,c);out(:,:)=NaN;out2=out;out3=out2;outh=out3;outv=outh;h = waitbar(0,'Please wait...');cpt=0;
for i=2:n-1;%boucle sur les latitudes
    for j=2:c-1;cpt=cpt+1;waitbar(cpt/((n-2).*(c-2)))%boucle sur les longitudes
        if isnan(mres(i,j))==0;
           Pa=mres(i+1,j-1);Pb=mres(i+1,j);Pj=mres(i+1,j+1);Pd=mres(i,j-1);Pf=mres(i,j+1);Pg=mres(i-1,j-1);Ph=mres(i-1,j);Pi=mres(i-1,j+1);
           Pa2=Pa;Pb2=Pb;Pj2=Pj;Pd2=Pd;Pf2=Pf;Pg2=Pg;Ph2=Ph;Pi2=Pi;
           if (isnan(Pa)==1);Pa2=0;end;if (isnan(Pb)==1);Pb2=0;end;if (isnan(Pj)==1);Pj2=0;end;if (isnan(Pd)==1);Pd2=0;end;
           if (isnan(Pf)==1);Pf2=0;end;if (isnan(Pg)==1);Pg2=0;end;if (isnan(Ph)==1);Ph2=0;end;if (isnan(Pi)==1);Pi2=0;end;
           ncell=[Pa,Pj,Pd,Pd,Pf,Pf,Pg,Pi];tfncell=~isnan(ncell);scell=sum(tfncell);%on cherche combien on a de cellule vide pour pouvoir les exclures
           Shor=((Pj2+(Pf2.*2)+Pi2)-(Pa2+(Pd2.*2)+Pg2))/(scell.*szlon);clear ncell tfncell scell;%Shor=((Pj+(Pf.*2)+Pi)-(Pa+(Pd.*2)+Pg))/(8.*szlon);
           %%modification de szlat ici en fonction de la taille des cellule due au projections!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           szlat=deg2km(1.0056,'earth');
           vlat=[wrdmulti(i+1,j-1),wrdmulti(i+1,j),wrdmulti(i+1,j+1),wrdmulti(i,j-1),wrdmulti(i,j+1),wrdmulti(i-1,j-1),wrdmulti(i-1,j),wrdmulti(i-1,j+1)];
                
           ncell=[Pa,Pb,Pb,Pj,Pg,Ph,Ph,Pi];tfncell=~isnan(ncell);scell=vlat(tfncell==1);
           Svert=((Pa2+(Pb2.*2)+Pj2)-(Pg2+(Ph2.*2)+Pi2))/(sum(scell.*szlat));clear ncell tfncell scell;%Svert=((Pa+(Pb.*2)+Pj)-(Pg+(Ph.*2)+Pi))/(sum(vlat.*szlat));
           %%======================================
           out(i,j)=sqrt((Shor.^2)+(Svert.^2));outh(i,j)=Shor;outv(i,j)=Svert;  
           out2(i,j)=atan(out(i,j)) * 57.29578;%on calcule l'heterogeneite spatiale
           out3(i,j)=atan2(Shor,Svert) * 57.29578;%on calcule le vecteur direction spatiale
           clear Pa Pb Pj Pd Pf Pg Ph Pi Shor Svert vlat;
        end;
    end;
end;close(h);[row,col]=find(out3<0);for k=1:length(row);out3(row(k),col(k))=360+out3(row(k),col(k));end;
save -ascii Spgrad_MegaExtinction(raster1x1)_v3.txt out2;save -ascii Bearing_MegaExtinction(raster1x1)_v3.txt out3;  
save -ascii Latgrad_MegaExtinction(raster1x1)_v3.txt outv;save -ascii Longrad_MegaExtinction(raster1x1)_v3.txt outh;

clat=flipud(yg');clon=xg;cutlat=find(clat==-26.5);cutlon=find(clon==134.5);out3(1:cutlat,:)=NaN;out3(:,1:cutlon)=NaN;out3(:,340:end)=NaN;
TF=isnan(out3);[row,col]=find(TF==0);nm=length(row);res2=NaN(nm,5);
for i=1:nm;res2(i,:)=[row(i),col(i),clat(row(i)),clon(col(i)),out3(row(i),col(i))];end;save -ascii Megafauna_SEBearing(LatLon)v3.txt res2;

clr=brewermap(12,'Paired');
figure(1),polarhistogram(res2(:,5),10,'Normalization','pdf','FaceColor',clr(7,:),'EdgeColor',clr(8,:));
          set(gca,'ThetaDir','clockwise','ThetaZeroLocation','top','LineWidth',2,'FontSize',13),rticklabels({});
          tr1=title('Megafauna Extinction Timning Bearing');set(tr1,'Fontsize',13);







  