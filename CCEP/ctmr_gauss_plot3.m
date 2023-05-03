function [elec_surf] = ctmr_gauss_plot3(cortex,electrodes,weights,cmap, hemi,viewside, do_lighting, gsp, clim, interp, show_colorbar, covgAdj, use_elecSurf )
% function [electrodes]=ctmr_gauss_plot(cortex,electrodes,weights)
% projects electrode locations onto their cortical spots in the
% left hemisphere and plots about them using a gaussian kernel
% for only cortex use:
% ctmr_gauss_plot(cortex,[0 0 0],0)
% rel_dir=which('loc_plot');
% rel_dir((length(rel_dir)-10):length(rel_dir))=[];
% addpath(rel_dir)

%     Copyright (C) 2009  K.J. Miller & D. Hermes, Dept of Neurology and Neurosurgery, University Medical Center Utrecht
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

%   Version 1.1.0, released 26-11-2009

% this codes is adapted from ctmr_gauss_plot for showing overlayed heatmap
% of two sets of electrodes (dian lyu 2022 Jul 14)
% weights is a matrix with more than one columns of numbers; requires
% covgAdj is whether to adjust for electrode coverage
% if covAdj is 1/2, it means that the hemispheres are grouped to one side,
% then the covAdj is only adjusted to the half.
% different colormaps

if nargin < 4 || isempty(cmap)
    cmap = [cbrewer2('seq','OrRd',64, 'spline')];%...
    %    cbrewer2('seq','RdPu',64, 'spline')]; % 

end

if isempty(weights)
    weights = ones(length(electrodes(:,1)),1);
end

if weights == 0
    weights = zeros(length(electrodes(:,1)),1);
end

m = size(cmap,1)/size(weights,2);

if nargin<5
    hemi = 'lh';
end
if nargin<6
    viewside = 'medial';
end

if nargin<7
    do_lighting=1;
end

if nargin < 8
    gsp = 50; % gaussian speading parameter
end

if nargin < 9
    clim = [];
end

if nargin < 10
    interp = 1; % methods of interpolation of the spatial spreading
end

if nargin < 11
    show_colorbar = 0;
end

if nargin < 12
    covgAdj = 0;

end

if nargin < 13
    use_elecSurf = 0; % whether use interpolated electrodes that are projected to the closest surface point
end
%load in colormap
%load('loc_colormap')
%load('loc_colormap_thresh')
% load('BlWhRdYl_colormap')
% load('BlGyOrCp_colormap')


brain=cortex.vert;

% %view from which side?
% temp=1;
% while temp==1
%     disp('---------------------------------------')
%     disp('to view from right press ''r''')
%     disp('to view from left press ''l''');
%     v=input('','s');
%     if v=='l'
%         temp=0;
%     elseif v=='r'
%         temp=0;
%     else
%         disp('you didn''t press r, or l try again (is caps on?)')
%     end
% end

if size(weights,1)~=length(electrodes(:,1))
    error('you sent a different number of weights than electrodes (perhaps a whole matrix instead of vector)')

end
%gaussian "cortical" spreading parameter - in mm, so if set at 10, its 1 cm
%- distance between adjacent electrodes

x = cortex.vert(:,1);
y = cortex.vert(:,2);
z = cortex.vert(:,3);
trids = cortex.tri;

%% plot
ax=[];
ax = newplot(ax);

axis tight
axis equal
hold on
if version('-release')>=12
    cameratoolbar('setmode', 'orbit')
else
    rotate3d on
end% if

%% find the closest points on the surface (DL)
elec_surf = nan(size(electrodes));

for i=1:size(electrodes,1)
    b_z=abs(brain(:,3)-electrodes(i,3));
    b_y=abs(brain(:,2)-electrodes(i,2));
    b_x=abs(brain(:,1)-electrodes(i,1));

    d = b_x.^2+b_z.^2+b_y.^2;
    elec_surf(i,:) = brain(d == min(d),:);
end

if  use_elecSurf == 1
    electrodes = elec_surf;
end

%% calculated gradient weights projected to the brain surface
for iw = 1:size(weights,2)
    %iw=1;
    weight = weights(:,iw);%+ m.*(iw-1);
    c = zeros(length(cortex(:,1)),1);
    c0 = zeros(length(cortex(:,1)),1);
    %cm = cmap{iw};
    for i=1:length(electrodes(:,1))
        b_z=abs(brain(:,3)-electrodes(i,3));
        b_y=abs(brain(:,2)-electrodes(i,2));
        b_x=abs(brain(:,1)-electrodes(i,1));
        if interp == 2
            d = exp((-(b_x.^2+b_z.^2+b_y.^2).^.5)/gsp^.5); %exponential fall off
            
        elseif interp == 1
            d = exp((-(b_x.^2+b_z.^2+b_y.^2))/gsp); 
        end
        % pure coverage
        c0 = c0+d';
        % weighted sum for activation
        d = weight(i).*d;
        c = c+d';
    end

    if covgAdj ~=0
%       col = c'./(c0'+1); % adjust for the elec coverage
      col =  c'./(sqrt(c0'.^2 + 1) .* covgAdj) ;
    
    else
        col = c';
    end
      

    if [1 3] == sort(size(col))
        col = repmat(col(:)', [size(cortex.vert, 1) 1]);
    end% if

%     if iw > 1
%         col(col<0.0001) = 0;
%     end
    
    h(iw) = patch('faces', trids,'vertices',[x(:) y(:) z(:)],...
        'FaceAlpha', 0.45, 'FaceColor', 'interp',  'FaceVertexCData', col,...
        'EdgeColor', 'none', 'EdgeAlpha', 0.0);%'parent',ax,
hold on
    %% find index to colormap
    % linear mapping
    cd = ones(size(col));
    [col_sorted, ind] = sort(col);

    if isempty(clim)
        cmax = max(col_sorted, [], 'omitnan');
        cmin = min(col_sorted, [], 'omitnan');
    else
        cmax = clim(2);
        cmin = clim(1);
    end
       col_rescaled = (col_sorted - cmin).*(m-1)/(cmax-cmin);
  
    cd(ind) = col_rescaled;
   
    
    cd = min(m.*iw, cd + m.*(iw-1));
    %cd = min(m,round(m*(col-min(col))/(max(col)-min(col)))+1) + m.*(iw-1);

    set(h(iw), 'CData',cd);
    
    shading interp;
%     a=get(gca);
%     %%NOTE: MAY WANT TO MAKE AXIS THE SAME MAGNITUDE ACROSS ALL COMPONENTS TO REFLECT
%     %%RELEVANCE OF CHANNEL FOR COMPARISON's ACROSS CORTICES
%     d = a.CLim;
%     set(gca,'CLim',[0 max(abs(d))])

%lighting phong; %play with lighting...
%material shiny;
material dull;
end % for iw

if isempty(clim)
    caxis([0, max(cd)])
else
    caxis([1, (iw.*m)])
end

colormap(gca, cmap((1:m) + m.*(iw-1), :))

if show_colorbar == 1

cbar = colorbar(gca);
cbar.Ticks = [1, 1+0.25.*(m-1),  1+0.5.*(m-1), 1+0.75.*(m-1), m];
cbar.TickLabels = {sprintf('%.2f', cmin), sprintf('%.2f', cmin + (cmax-cmin).*0.25), ...
    sprintf('%.2f', cmin + (cmax-cmin).*0.5), sprintf('%.2f', cmin + (cmax-cmin).*0.75), ...
    sprintf('%.2f', cmax)};
cbar.FontSize = 12;
%cbar.Position = cbar.Position ;

end

%end
%colormap(jet);
% view(ax,3);

%material([.3 .8 .1 10 1]);
%     material([.2 .9 .2 50 1]); %  BF: editing mesh viewing attributes
switch ax.NextPlot
    case {'replaceall','replace'}
        view(ax,3);
        grid(ax,'on');
    case {'replacechildren'}
        view(ax,3);
end

ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

axis off
%%
%set(gcf,'Renderer', 'zbuffer','Position',[500 500 900 900]); % BF: added for lateral. view
%set(gcf,'Renderer', 'zbuffer','Position',[400 400 500 900]); % BF: added for inf. view
%set(gcf,'Renderer', 'zbuffer','Position',[400 400 950 550]); % BF: added figure size for movie


%% decorations

% if v=='l'
if do_lighting
    l=light;

    if strcmp(hemi,'lh')
        view(270, 0);
        % set(l,'Position',[-1 0 1])
        set(l,'Position',[-1 0 0],'Color',[0.8 0.8 0.8]);
        % elseif v=='r'
    elseif strcmp(hemi,'rh')
        view(90, 0);
        % set(l,'Position',[1 0 1])
        set(l,'Position',[1 0 0],'Color',[0.8 0.8 0.8]);
    end
end


if strcmp(hemi,'left')
    switch viewside
        case 'medial'
            loc_view(90, 0)
        case 'lateral'
            loc_view(270, 25)
        case 'anterior'
            loc_view(180,0)
        case 'posterior'
            loc_view(0,0)
        case 'ventral'
            loc_view(180,270)
        case 'temporal'
            loc_view(270,-80)
        case 'dorsal'
            loc_view(0,90)
        case 'latero-ventral'
            loc_view(270,-45)
        case 'medio-dorsal'
            loc_view(90,45)
        case 'medio-ventral'
            loc_view(90,-45)
        case 'medio-posterior'
            loc_view(45,0)
        case 'medio-anterior'
            loc_view(135,0)
        case 'frontal'
            loc_view(-120,10)
        case 'parietal'
            loc_view(-70,10)
    end
    %     set(l,'Position',[-1 0 1])
elseif strcmp(hemi,'right')
    switch viewside
        case 'medial'
            loc_view(272, 0)
        case 'lateral'
            loc_view(90, 25)
        case 'anterior'
            loc_view(180,0)
        case 'posterior'
            loc_view(0,0)
        case 'dorsal'
            loc_view(0,90)
        case 'ventral'
            loc_view(180,270)
        case 'temporal'
            loc_view(90,-80)
        case 'latero-ventral'
            loc_view(90,-45)
        case 'medio-dorsal'
            loc_view(270,45)
        case 'medio-ventral'
            loc_view(270,-45)
        case 'medio-posterior'
            loc_view(315,0)
        case 'medio-anterior'
            loc_view(225,0)
        case 'frontal'
            loc_view(120,10)
        case 'parietal'
            loc_view(70,10)
    end

else
    error('hemisphere should be either left or right')
    %     set(l,'Position',[1 0 1])
end

set(gcf,'color','w')

% hcb = colorbar;
% set(hcb,'YTick',[])

% view(270, 0);
% set(l,'Position',[-1 0 1])
% elseif v=='r'
% view(90, 0);
% set(l,'Position',[1 0 1])
% end
% %exportfig
% exportfig(gcf, strcat(cd,'\figout.png'), 'format', 'png', 'Renderer', 'painters', 'Color', 'cmyk', 'Resolution', 600, 'Width', 4, 'Height', 3);
% disp('figure saved as "figout"');
axis off

% %exportfig
% exportfig(gcf, strcat(cd,'\figout.png'), 'format', 'png', 'Renderer', 'painters', 'Color', 'cmyk', 'Resolution', 600, 'Width', 4, 'Height', 3);
% disp('figure saved as "figout"');

