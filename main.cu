#include "toolpath.h"

void handle_error_cuda(char *line){
    cudaError_t err = cudaGetLastError();
    if( cudaSuccess != err){
        printf("%s) Cuda error: %s.\n", line, cudaGetErrorString( err) );
        exit(0);
    }                        
}

__device__ double dMin (double a, double b) {
	double output = a;
	if (a>b) output = b;
	if (b>a) output = a;
	return output;
}

__device__ double dMax (double a, double b) {
	double output = a;
	if (a>b) output = a;
	if (b>a) output = b;
	return output;
}
__global__ void checkIntersection(double *lineX0, double *lineY0, double *lineX1, double *lineY1, short *isThereIntersection, long int numOfLines, long int lineId, double *catchedIntersectionX, double *catchedIntersectionY, double eps){
	int idx = blockDim.x* blockIdx.x + threadIdx.x;
	if (idx < numOfLines) {
		if (idx != lineId) {
			isThereIntersection[idx] = 0;
			if ((lineX0[idx] != lineX1[idx]) && (lineX0[lineId] != lineX1[lineId])) {
				double line1_k = (lineY1[idx] - lineY0[idx]) 		/ (lineX1[idx] - lineX0[idx]);
				double line2_k = (lineY1[lineId] - lineY0[lineId]) 	/ (lineX1[lineId] - lineX0[lineId]);

				double line1_y0 = lineY0[idx] 		- line1_k * lineX0[idx];
				double line2_y0 = lineY0[lineId] 	- line2_k * lineX0[lineId];
				if (line2_k != line1_k) {
					double intersection_x = (line1_y0 - line2_y0) / (line2_k - line1_k);
					double intersection_y = line1_y0 + line1_k * intersection_x;
					if ((intersection_x >= dMin(lineX0[idx], 	lineX1[idx]) - eps) 			&&
						(intersection_x <= dMax(lineX0[idx], 	lineX1[idx]) + eps) 			&&
						(intersection_x >= dMin(lineX0[lineId],  lineX1[lineId]) - eps) 		&&	
						(intersection_x <= dMax(lineX0[lineId], 	lineX1[lineId]) + eps)	&&
						(intersection_y >= dMin(lineY0[idx], lineY1[idx]) - eps)				&&
						(intersection_y <= dMax(lineY0[idx], lineY1[idx]) + eps)				&&
						(intersection_y >= dMin(lineY0[lineId], lineY1[lineId]) - eps)		&&
						(intersection_y <= dMax(lineY0[lineId], lineY1[lineId]) + eps)		
						)
					   
					{
						isThereIntersection[idx] = 1;
						catchedIntersectionX[idx] = intersection_x;
						catchedIntersectionY[idx] = intersection_y;
					}
				}
			}
		}
	}
}
#include "detectOuter.h"

short  sign_host(double x) {
	short output = 0;
	if (x<0) output = -1;
	if (x>0) output = 1;
	return output;
}

double vetctorMult_host(double ax, double ay, double bx, double by) {
	return ax * by - ay * bx;
}

int checkPoint_host(double dotX, double dotY, double *polygonX, double *polygonY, long int *profileLength, long int *profileStartIndex, long int numOfProfiles, long int profileId){
	long int startIndex = profileStartIndex[profileId];
	double abX = polygonX[startIndex+1] - polygonX[startIndex];
	double abY = polygonY[startIndex+1] - polygonY[startIndex];
	double apX = dotX - polygonX[startIndex];
	double apY = dotY - polygonY[startIndex];
	double bpX = dotX - polygonX[startIndex+1];
	double bpY = dotY - polygonY[startIndex+1];
	short pointIn = 1;
	short currSign = sign_host(vetctorMult_host(abX,abY,apX,apY));
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
			currSign = sign_host(vetctorMult_host(abX,abY,apX,apY));
		} else {
			abX = polygonX[startIndex] - polygonX[i];
			abY = polygonY[startIndex] - polygonY[i];
			apX = dotX - polygonX[i];
			apY = dotY - polygonY[i];
			bpX = dotX - polygonX[startIndex];
			bpY = dotY - polygonY[startIndex];
			currSign = sign_host(vetctorMult_host(abX,abY,apX,apY));
		}
		
		if (currSign * prevSign < 0) {pointIn = 0; break;}
	}
	return pointIn;
}

int main(int argc, char const *argv[]){
	#include "inputData"
	#include "detectLines"
	#include "detectProfiles"
	#include "copyDataToDevice"
	#include "printLines"
	#include "detectIntersections"
	// #include "printIntersections"
	#include "copyIntersectionsToCUDA"
	// for (int i = 0; i < numOfPoints; i++){
	// 	printf("%lf %lf %ld\n",polygonX_host[i],polygonY_host[i],polygonId_host[i] );
	// }
	// double *finalIntersectionsX = (double *)malloc(sizeof(double));
	// double *finalIntersectionsY = (double *)malloc(sizeof(double));
	// long int tempNumOfIntersections = 0;
	// for (long int i = 0; i < numOfIntersections; i++) {
	// 	int numOfRepeates = 0;
	// 	for (int j = i+1; j < numOfIntersections; j++){
	// 		if ((fabs(intersectionsX_host[i] - intersectionsX_host[j]) <= eps) &&
	// 			(fabs(intersectionsY_host[i] - intersectionsY_host[j]) <= eps)
	// 			)
	// 		numOfRepeates++;
	// 		intersectionsY_host[j] = 
	// 	}
	// 	if ()
	// }
	#include "detectOuter"
	
	
	
	FILE *outer = fopen("intersections.txt","w");
	for (long int i = 0; i < numOfIntersections; i++) {
		// pointIsOuter_host[i]=1;
		if (pointIsOuter_host[i] >= 1) {
			numOfEdgePoints++;
			fprintf(outer, "%lf,%lf\n",intersectionsX_host[i],intersectionsY_host[i] );
		}
	}
	// float radius = 5;
	// for (float angle = 0; angle < 360; angle+=0.001){
	// 	float tempX = radius * cos(angle * 3.14 / 180);
	// 	float tempY = radius * sin(angle * 3.14 / 180);
	// 	tempX = (float)(random() % 2000)/100 - 8;
	// 	tempY = (float)(random() % 5000)/100 - 20;
	// 	for (int i = 0; i < numOfProfiles; i++) {
	// 		if (checkPoint_host(tempX, tempY, polygonX_host, polygonY_host, profileLength_host, profileStartIndex_host, numOfProfiles, i) == 1) {
	// 			int numOfIntersectedProfiles = 0;
	// 			float rad = 0.1;
	// 			for (float x = 0; x < 360; x+=20) {
	// 				double tx = tempX + rad * cos(x*3.14/180);
	// 				double ty = tempY + rad * sin(x*3.14/180);
	// 				if (checkPoint_host(tx, ty, polygonX_host, polygonY_host, profileLength_host, profileStartIndex_host, numOfProfiles, i) == 1) {
	// 					numOfIntersectedProfiles++;
	// 				}
	// 			}
	// 			if (numOfIntersectedProfiles<=1) fprintf(outer, "%lf,%lf\n",tempX,tempY );
	// 		}
	// 	}
	// }
	
	fclose(outer);

	cudaFree(polygonX_device);
	cudaFree(polygonY_device);
	cudaFree(polygonId_device);
	cudaFree(lineX0_device);
	cudaFree(lineY0_device);
	cudaFree(lineX1_device);
	cudaFree(lineY1_device);
	cudaFree(catchedIntersectionX_device);
	cudaFree(catchedIntersectionY_device);
	cudaFree(isThereIntersection_device);
	cudaFree(intersectionsX_device);
	cudaFree(intersectionsY_device);
	cudaFree(profileStartIndex_device);
	cudaFree(profileLength_device);
	cudaFree(pointIsOuter_device);

	free(pointIsOuter_host);
	free(profileLength_host);
	free(profileStartIndex_host);
	free(isThereIntersection_host);
	free(lineX0_host);
	free(lineY0_host);
	free(lineX1_host);
	free(lineY1_host);
	free(catchedIntersectionX_host);
	free(catchedIntersectionY_host);
	free(polygonX_host );
	free(polygonY_host );
	free(polygonId_host);
	free(intersectionsX_host);
	free(intersectionsY_host);


	printf("numOfPoints=%ld\n",numOfPoints );
	printf("numOfLines=%ld\n",numOfLines );
	printf("numOfProfiles=%ld\n",numOfProfiles );
	printf("numOfIntersections=%ld\n",numOfIntersections );
	printf("numOfEdgePoints=%ld\n", numOfEdgePoints );
	return 0;
}