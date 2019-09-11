function txt = myupdatefcn(~,event_obj,ID,index)
% Customizes text of data tips
pos = get(event_obj,'Position');
I = get(event_obj, 'DataIndex');
txt = {['X: ',num2str(pos(1))],...
    ['Y: ',num2str(pos(2))],...     
    ['i: ',num2str(index(I))],...
    ['ID: ',num2str(ID(I))]};

end