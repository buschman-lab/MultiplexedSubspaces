function y = movmode(x,w)

%perform a sliding window mode, centered at x +/-w;

y = NaN(size(x));
x = cat(1,NaN(w,1),x,NaN(w,1));

window = (-w:w);
COUNT=1;
for i = (w+1):(numel(x)-w)
    y(COUNT) = mode(x(i+window));
    COUNT = COUNT+1;
end
   