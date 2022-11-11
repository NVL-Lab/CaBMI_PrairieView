function cp = getCurrentPoint(hAx)
    cp = get(hAx, 'CurrentPoint');
    cp = cp(1, 1:2);
end