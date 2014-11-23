__device__ short  sign_device(double x) {
	short output = 0;
	if (x<0) output = -1;
	if (x>0) output = 1;
	return output;
}

__device__ double vetctorMult_device(double ax, double ay, double bx, double by) {
	return ax * by - ay * bx;
}

__device__ int checkPoint_device(double dotX, double dotY, double *polygonX, double *polygonY, long int *profileLength, long int *profileStartIndex, long int numOfProfiles, long int profileId){
	long int startIndex = profileStartIndex[profileId];
	double abX = polygonX[startIndex+1] - polygonX[startIndex];
	double abY = polygonY[startIndex+1] - polygonY[startIndex];
	double apX = dotX - polygonX[startIndex];
	double apY = dotY - polygonY[startIndex];
	double bpX = dotX - polygonX[startIndex+1];
	double bpY = dotY - polygonY[startIndex+1];
	short pointIn = 1;
	short currSign = sign_device(vetctorMult_device(abX,abY,apX,apY));
	for (long int i = startIndex + 1; i < startIndex + profileLength[profileId] ; i++){
		short prevSign = currSign;
		if (prevSign == 0) {
			if 	(((bpX*bpX + bpY*bpY) > (abX*abX + abY * abY)) ||
				((apX*apX + apY*apY) > (abX*abX + abY * abY))) 
			{
				 pointIn = 0;
				 break;
			}
		}
		
		if (i != startIndex + profileLength[profileId] - 1 ){
			abX = polygonX[i+1] - polygonX[i];
			abY = polygonY[i+1] - polygonY[i];
			apX = dotX - polygonX[i];
			apY = dotY - polygonY[i];
			bpX = dotX - polygonX[i+1];
			bpY = dotY - polygonY[i+1];
			currSign = sign_device(vetctorMult_device(abX,abY,apX,apY));
		} else {
			abX = polygonX[startIndex] - polygonX[i];
			abY = polygonY[startIndex] - polygonY[i];
			apX = dotX - polygonX[i];
			apY = dotY - polygonY[i];
			bpX = dotX - polygonX[startIndex];
			bpY = dotY - polygonY[startIndex];
			currSign = sign_device(vetctorMult_device(abX,abY,apX,apY));
		}
		
		if (currSign * prevSign < 0) {pointIn = 0; break;}
	}
	return pointIn;
}

// __global__ void detectOuter(long int intersectionStartIndex, short *pointIsOuter, long int numOfIntersections, double *intersectionsX, double *intersectionsY,double *polygonX, double *polygonY, long int *profileLength, long int *profileStartIndex, long int numOfProfiles){
// 	long int idx = blockDim.x* blockIdx.x + threadIdx.x + intersectionStartIndex;
// 	// printf("idx=%ld\n", idx);
// 	if (idx < numOfIntersections) {
// 		pointIsOuter[idx] = 1;
// 		double intersectionX = intersectionsX[idx];
// 		double intersectionY = intersectionsY[idx];
// 		double radius = 0.1;
// 		for (float angle = 0; angle < 360; angle+=30) {
// 			int numOfIntersectedProfiles = 0;
// 			double tempX = intersectionX + radius * cos(angle * 3.14 / 180);
// 			double tempY = intersectionY + radius * sin(angle * 3.14 / 180);
// 			for (long int polygonIndex = 0; polygonIndex < numOfProfiles; polygonIndex++) {
// 				if (checkPoint_device(tempX, tempY, polygonX, polygonY, profileLength, profileStartIndex, numOfProfiles, polygonIndex) == 1) {
// 					numOfIntersectedProfiles++;
// 					break;
// 				}
// 			}
// 			if (numOfIntersectedProfiles<1) pointIsOuter[idx] = 1;
// 		}
// 	}
// }

// __global__ void checkProfiles_device(short *pointIsOuter, long int pointId, double x, double y, double *polygonX, double *polygonY, long int numOfProfiles, long int *profileLength, long int *profileStartIndex) {
// 	long int idx = blockDim.x* blockIdx.x + threadIdx.x;
// 	if (idx < numOfProfiles) {
// 		if (pointIsOuter[pointId] == 0 ) {
// 			int numOfIntersectedProfiles = 0;
// 			for (long int polygonIndex = 0; polygonIndex < numOfProfiles; polygonIndex++) {
// 				if (checkPoint_device(x, y, polygonX, polygonY, profileLength, profileStartIndex, numOfProfiles, polygonIndex) == 1) {
// 					numOfIntersectedProfiles++;
// 					break;
// 				}
// 			}
// 			if (numOfIntersectedProfiles<1) pointIsOuter[pointId]++;
// 		}
// 	}
// }

__global__ void checkProfiles_device(short *pointIsOuter, long int pointId, double x, double y, double *polygonX, double *polygonY, long int numOfProfiles, long int *profileLength, long int *profileStartIndex) {
	long int idx = blockDim.x* blockIdx.x + threadIdx.x;
	if (idx < numOfProfiles) {
		if (pointIsOuter[pointId] == 0 ) {
			float radius = 0.1;
			for (float angle = 0; angle < 360; angle+=45) {
				int numOfIntersectedProfiles = 0;
				double tX = x + radius * cos(angle * 3.14 / 180);
				double tY = y + radius * sin(angle * 3.14 / 180);
				for (long int polygonIndex = 0; polygonIndex < numOfProfiles; polygonIndex++) {
					if (checkPoint_device(tX, tY, polygonX, polygonY, profileLength, profileStartIndex, numOfProfiles, polygonIndex) == 1) {
						numOfIntersectedProfiles++;
						break;
					}
				}
				if (numOfIntersectedProfiles<1) {
					pointIsOuter[pointId]++; 
					break;
					asm("exit;");
				}
			}
			
		}
	}
}
