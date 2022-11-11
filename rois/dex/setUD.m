function out = setUD(handle, paramname, value)
    u = get(handle, 'UserData');
    u.(paramname) = value;
    set(handle, 'UserData', u);