function [] =setXlabelsImage(data,xticklabels)
                xticks = linspace(1, size(data, 2), numel(xticklabels));
                set(gca, 'XTick', xticks, 'XTickLabel', round(xticklabels,2))
