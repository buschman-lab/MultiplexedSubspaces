function uC = UniqueCell(C)
% Reply cell array with unique elements, independent of the type
% and size of the cell elements.
% Author: Jan Simon, Heidelberg, License: CC BY-SA 3.0
u = true(size(C));
for iC = 1:numel(C)
  if u(iC)
    aC = C{iC};
    for jC = iC + 1:numel(C)
      u(jC) = u(jC) && ~isequal(aC, C{jC});
    end
  end
end
uC = C(u);
end