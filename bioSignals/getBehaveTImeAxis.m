function behave_events_times = getBehaveTImeAxis(t, bejaveLabels, del)
if ~isempty(intersect(unique(bejaveLabels), 9))
    
    s = hist(bejaveLabels',0:9);
m = mean(s,2);
st = std(s, 1, 2);
th = m+del*st;
globalbehavelabel = zeros(length(t),1);
inds = find(s([ 3],:)>th(3));
while ~all(diff(inds)==1)
inds = inds((diff(inds)==1));
end
globalbehavelabel(inds) = 3;

inds = find(s([ 2],:)>th(2));
while ~all(diff(inds)==1)
inds = inds((diff(inds)==1));
end
globalbehavelabel(inds) = 2;

inds = find(s([ 5],:)>th(5));
while ~all(diff(inds)==1)
inds = inds((diff(inds)==1));
end
globalbehavelabel(inds) = 4;

inds = find(s([ 10],:)>th(10));
while ~all(diff(inds)==1)
inds = inds((diff(inds)==1));
end
globalbehavelabel(inds) = 9;

behave_events_times = t((diff(globalbehavelabel)~=0));
else
s = hist(bejaveLabels',0:4);
m = mean(s,2);
st = std(s, 1, 2);
th = m+del*st;
globalbehavelabel = zeros(length(t),1);
inds = find(s([ 3],:)>th(3));
while ~all(diff(inds)==1)
inds = inds((diff(inds)==1));
end
globalbehavelabel(inds) = 3;

inds = find(s([ 2],:)>th(2));
while ~all(diff(inds)==1)
inds = inds((diff(inds)==1));
end
globalbehavelabel(inds) = 2;

inds = find(s([ 5],:)>th(5));
while ~all(diff(inds)==1)
inds = inds((diff(inds)==1));
end
globalbehavelabel(inds) = 4;


behave_events_times = t((diff(globalbehavelabel)~=0));
end