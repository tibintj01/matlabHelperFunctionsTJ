function [] = addCustomLegendTibin(colors,names)
	%axes( 'Position', [0.95, 0.95, 0.05, 0.05] )
axes( 'Position', [0.99, 0.99, 0.01, 0.01] )

	numEntries=length(colors);
	for i=1:numEntries
        try
            plot(NaN,NaN,[colors{i}]);
        %catch
        %    plot(NaN,NaN,[colors(i,:)]);
        end
		hold on
	end

	legend(names{:})

        ax = gca
        ax.Visible = 'off'
