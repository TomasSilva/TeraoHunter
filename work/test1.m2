needsPackage "HyperplaneArrangements"
R = QQ[x,y,z,w];
A = arrangement {x,y,z,w,x-y,x-z,x-w,y-z,y-w,z-w};
D = der A;
print toString betti res prune image D;
print isFreeModule prune image D;
print toString degrees source D;
