function make_network_bargraph(networkIDs,networkvals,errors,sort_byval)
%make_network_bargraph(networkIDs,networkvals,[errors],[sort_byval])

if exist('errors') && isempty(errors)
    clear errors
end

if ~exist('sort_byval')
    sort_byval = true;
end


networkIDs(isnan(networkvals)) = [];
networkvals(isnan(networkvals)) = [];
networkIDs(networkIDs<1) = 18;


power_surf_colormap = [1 0 0;0 0 .8;1 1 0;1 .8 .6;0 1 0;1 .6 1;0 .6 .6;0 0 0;.35 0 .65;.2 1 1;1 .5 0;.65 .25 1;0 .25 .6;.6 1 .6;.2 .3 1;1 1 1;0 .4 0; repmat([.25 .25 .25],50,1)];

[~,sorti] = sort(networkvals,'descend');
if ~sort_byval
    sorti = 1:length(networkvals);
end


figure;
set(gcf,'Position',[813 30 1102 805])
set(gcf,'Color',[1 1 1]);%set(gcf,'Color',[.9 .9 .9])
set(gca,'Color',[1 1 1]);%set(gca,'Color',[.9 .9 .9])

hold on

for networknum = 1:length(networkIDs);
    h = bar(find(sorti==networknum),networkvals(networknum));
    
    decimalval = mod(networkIDs(networknum),1);
    if decimalval==0
        thiscolor = power_surf_colormap(networkIDs(networknum),:);
    else
        thiscolor = sum([power_surf_colormap(floor(networkIDs(networknum)),:) .* (1-decimalval) ; power_surf_colormap(ceil(networkIDs(networknum)),:) .* (decimalval)],1);
    end
    set(h,'FaceColor', thiscolor);
    
    if exist('errors','var')
        he = errorbar(find(sorti==networknum),networkvals(networknum),errors(networknum),'k');
        
        if verLessThan('matlab','8.5')
            
            set(he,'LineStyle','none');
            hb = get(he,'children');
            Xdata = get(hb(2),'Xdata');
            temp = 4:3:length(Xdata);
            temp(3:3:end) = [];
            % xleft and xright contain the indices of the left and right
            %  endpoints of the horizontal lines
            xleft = temp; xright = temp+1;
            % Increase line length by 0.2 units
            Xdata(xleft) = Xdata(1) - .1;
            Xdata(xright) = Xdata(1) + .1;
            set(hb(2),'Xdata',Xdata)
            
        end
        
    end
end

% if exist('errors','var')
%         he = errorbar(sorti,networkvals,errors,'k');
%         set(he,'LineStyle','none');
%         end
    
        
%         he = errorbar(find(sorti==networknum),networkvals(networknum),errors(networknum),'-k');
%         hesub = get(he,'children');
%         errorXs = get(hesub(2),'XData');
%         errorXs([4 7]) = errorXs(1)-.1;
%         errorXs([5 8]) = errorXs(1)+.1;
%         set(hesub(2),'XData',errorXs);



set(gca,'FontSize',30);
xlim([0 (length(networkIDs)+.6)])
set(gca,'XTick',[])
