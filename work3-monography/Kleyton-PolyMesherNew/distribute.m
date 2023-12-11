function df = distribute(P,Wells)
	NElem = 100;
    NWells = size(Wells,1);
    wx = Wells(:,1);
	wy = Wells(:,2);
    df = 0;
	for i = 1:NWells
		df = df + exp(-(NElem/100)*((((P(:,1)-wx(i))).^2)+((P(:,2)-wy(i)).^2)));
	end
end