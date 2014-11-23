struct point 		*getPoints(char const *filename) {
	FILE *file 			= fopen(filename,"r");
	if (file == NULL) {printf("No file!\n"); exit(0);}
	struct point *data 	= (struct point *)malloc(sizeof(struct point));
	int currIndex 		= 0;
	int currProfileId 	= -1;
	do{
		int gMode=0;
		fscanf(file,"%ld G%dX%lfY%lf\n",&data[currIndex].profileId,&gMode,&data[currIndex].x,&data[currIndex].y);		//printf("%ld %lf %lf\n",data[currIndex].profileId,data[currIndex].x,data[currIndex].y );
		data[currIndex].outer=0;
		if (currProfileId != data[currIndex].profileId) {
			currProfileId=data[currIndex].profileId;
		}
		currIndex++;
		numOfPoints++;
		data = (struct point *)realloc(data,(currIndex+1)*sizeof(struct point));
	}while (!feof(file));
	fclose(file);
	return data;
}
