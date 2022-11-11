function out = getUD(handle, paramname)
    u = get(handle, 'UserData');
    if nargin < 2
        out = u;
    else
        out = u.(paramname);
    end