

%s1=plot(t,dSlow,'k','LineWidth',3)
%hold on

dPerField=fieldCentersPerField;
tPerField=NaN(size(dPerField));
numFields=length(dPerField);

if(useHighSpeed)
    s2=plot(t,dFast,'k--','LineWidth',3);
    for i=1:numFields
        tPerField(i)=interp1(dFast,t,dPerField(i));
    end
    lineTypeStr='--';
    markerTypeStr='o';
      markerSize=15/2;
else
    s1=plot(t,dSlow,'k','LineWidth',3);
    for i=1:numFields
        tPerField(i)=interp1(dSlow,t,dPerField(i));
    end

    lineTypeStr='-';
     markerTypeStr='.';
     markerSize=50/2;
end


%h1=hline(d1,pastelGreen);
%h2=hline(d2,pastelPurple);

axis square
ylim([0 yMaxDisp])
xlim([0 xMaxDisp])
xlabel('Time')
ylabel('Distance')


if(useLogScaleTime)
     title('Distance vs Time (log-log scale)')
     set(gca,'xscale','log')
      set(gca,'yscale','log')
else
    title('Distance vs Time')
end

hold on
%{
plot([0.0001 t1],[d1 d1],lineTypeStr,'Color',pastelGreen,'LineWidth',3)
hold on
plot([0.0001 t2],[d2 d2],lineTypeStr,'Color',pastelPurple,'LineWidth',3)

plot([t1 t1],[0 d1],lineTypeStr,'Color',pastelGreen,'LineWidth',3)
hold on
plot([t2 t2],[0 d2],lineTypeStr,'Color',pastelPurple,'LineWidth',3)
%}

for i=1:numFields
    plot(tPerField(i),dPerField(i),markerTypeStr,'Color',fieldColors(i,:),'LineWidth',3,'MarkerSize',markerSize)
    hold on
    plot([0.0001 xMaxDisp],[dPerField(i) dPerField(i)],lineTypeStr,'Color',fieldColors(i,:),'LineWidth',3)
end

%{
plot(t1,d1,markerTypeStr,'Color',pastelGreen,'LineWidth',3,'MarkerSize',markerSize)
hold on
plot(t2,d2,markerTypeStr,'Color',pastelPurple,'LineWidth',3,'MarkerSize',markerSize)

plot([0.0001 xMaxDisp],[d1 d1],lineTypeStr,'Color',pastelGreen,'LineWidth',3)
hold on
plot([0.0001 xMaxDisp],[d2 d2],lineTypeStr,'Color',pastelPurple,'LineWidth',3)
%}

xlim([minLogTdisp maxLogTdisp])
ylim([minLogTdisp maxLogTdisp])
%ylim([0 1])
box off