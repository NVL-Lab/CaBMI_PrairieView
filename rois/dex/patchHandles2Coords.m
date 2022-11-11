function coords = patchHandles2Coords(patchHandles)

coords = cell(size(patchHandles));

for i = 1:numel(patchHandles)
    for j = 1:length(patchHandles{i})
        coords{i}(j) = {[get(patchHandles{i}(j), 'XData') get(patchHandles{i}(j), 'YData')]};
    end
end