function OK = isfigure(h)
%Checks if input is a figure handle or not;
%@MSeifert 2020


if strcmp(get(h,'type'),'figure')
  OK = true;
else
  OK = false;
end