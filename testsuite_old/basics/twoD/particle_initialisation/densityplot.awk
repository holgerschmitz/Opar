BEGIN{
    width = 200;
    height = 200;
    minx = 0;
    maxx = 1e-5;
    miny = 0;
    maxy = 1e-5;
    dx = (maxx-minx)/width;
    dy = (maxy-miny)/height;

    for (i=0; i<width; i++)
	for (j=0; j<height; j++)
	    den[i + width*j]=0;
}

{
    xi = int(($1-minx)/dx);
    yi = int(($2-miny)/dy);
    if ((xi>=0) && (xi<width) && (yi>=0) && (yi<height)) {
	den[xi + width*yi] += $3;
    }
}

END{
     for (j=0; j<height; j++) {
        for (i=0; i<width; i++)
	    printf("%g ",den[i + width*j]);
	printf("\n");
    }
}