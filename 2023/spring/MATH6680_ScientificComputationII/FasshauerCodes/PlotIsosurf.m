% PlotIsosurf(xe,ye,ze,Pf,neval,exact,maxerr,isomin,
%             isostep,isomax,caption)
% Generates plot of isosurfaces of Pf false colored by 
% the error abs(Pf-exact)
% isomin,isostep,isomax define the range and number of 
% isosurfaces, 
% and the caption is displayed.
function PlotIsosurf(xe,ye,ze,Pf,neval,exact,maxerr,...
         isomin,isostep,isomax,caption)
% Plot isosurfaces 
figure
hold on
for isovalue=isomin:isostep:isomax
    pfit = patch(isosurface(xe,ye,ze,reshape(Pf,neval,...
        neval,neval),isovalue,reshape(abs(Pf-exact),...
        neval,neval,neval)));
    isonormals(xe,ye,ze,reshape(Pf,neval,neval,neval),pfit)
    set(pfit,'FaceColor','interp','EdgeColor','none');
    daspect([1 1 1])
    view(3); axis([0 1 0 1 0 1])
end
[cmin cmax] = caxis;
caxis([cmin-.25*cmax cmax])
colormap hsv
vcb = colorbar('vert');
ylim(vcb,[0 cmax])
set(get(vcb,'YLabel'),'String','Error')
title(caption)
hold off

