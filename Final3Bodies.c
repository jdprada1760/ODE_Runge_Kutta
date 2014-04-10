#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#define PI 3.14159265358979323846
#define FLOAT double
#define solar_mass 1.0
#define G 39.486
#define T 5000.0
#define R 100.0
#define n_dim 2
#define n_masses 3
#define theta 2*PI/3.0
#define r R*sqrt(3.0)
#define v sqrt((1.1*G*solar_mass)/(3.0*r))

FLOAT norm( FLOAT vec1, FLOAT vec2 );
FLOAT dotP( FLOAT vec11, FLOAT vec12, FLOAT vec21, FLOAT vec22);
FLOAT * crossP( FLOAT vec1, FLOAT vec2);
FLOAT * allocate( int n );
void initPos( FLOAT *xX, FLOAT *xY );
void initVel( FLOAT *vX, FLOAT *vY );
void initMass( FLOAT *m );
void updateAcc( FLOAT *xX, FLOAT *xY, FLOAT *aX, FLOAT *aY , FLOAT *m);
void rungeKutta( FLOAT *xX, FLOAT *xY, FLOAT *vX, FLOAT *vY, FLOAT *m ,FLOAT *aX, FLOAT *aY, FLOAT dt);
FLOAT calcEnergy(FLOAT *x, FLOAT *y, FLOAT *vx, FLOAT *vy, FLOAT *m);

int main(int argc, char **argv){
 
  //------------------------------------------------------------------------
  // Constantes
  //------------------------------------------------------------------------

  // Archivo en formato t x y vx vy ax ay energia
  FILE *data1;
  data1 = fopen("3bodyProblem_1.data", "w");
  FILE *data2;
  data2 = fopen("3bodyProblem_2.data", "w");
  FILE *data3;
  data3 = fopen("3bodyProblem_3.data", "w");
  FILE *energydat;
  energydat = fopen("energy.data", "w");
  
  //------------------------------------------------------------------------
  // Variables 
  //------------------------------------------------------------------------
   
  // Crea Variables
  FLOAT *x;
  FLOAT *y;
  FLOAT *vX;
  FLOAT *vY;
  FLOAT *aX;
  FLOAT *aY;
  FLOAT *m;
  FLOAT actualTime = 0;
  FLOAT dt = atof(argv[1]);
  FLOAT energy = 0;
  // Guarda Memoria
  x = allocate(n_masses);
  y = allocate(n_masses);
  vX = allocate(n_masses);
  vY = allocate(n_masses);
  aX = allocate(n_masses);
  aY = allocate(n_masses);
  m = allocate(n_masses);
  // Inicializa Variables
  initPos(x, y);
  initVel(vX, vY);
  initMass(m);
  updateAcc(x,y,aX,aY,m);
  // Escribe en el archivo los primeros términos
  fprintf(data1, "%f %f %f %f %f %f %f\n",actualTime,x[0],y[0],vX[0],vY[0],aX[0],aY[0]);
  fprintf(data2, "%f %f %f %f %f %f %f\n",actualTime,x[1],y[1],vX[1],vY[1],aX[1],aY[1]);
  fprintf(data3, "%f %f %f %f %f %f %f\n",actualTime,x[2],y[2],vX[2],vY[2],aX[2],aY[2]);
  fprintf(energydat, "%f %f\n",actualTime,energy);
  
  do{

    //------------------------------------------------------------------------
    // 1st second order DE Runge - Kutta (x1,v1)
    //------------------------------------------------------------------------
    
    // Actualiza la aceleración (función de segunda derivada)
    updateAcc(x,y,aX,aY,m);
    // RungeKutta 4th order
    rungeKutta(x,y,vX,vY,m,aX,aY,dt);
    // Actualiza el tiempo
    actualTime += dt;
    // Calcula la energía del sistema
    energy = calcEnergy(x,y,vX,vY,m);
    // Escribe en el archivo
    fprintf(data1, "%f %f %f %f %f %f %f\n",actualTime,x[0],y[0],vX[0],vY[0],aX[0],aY[0]);
    fprintf(data2, "%f %f %f %f %f %f %f\n",actualTime,x[1],y[1],vX[1],vY[1],aX[1],aY[1]);
    fprintf(data3, "%f %f %f %f %f %f %f\n",actualTime,x[2],y[2],vX[2],vY[2],aX[2],aY[2]);
    fprintf(energydat, "%f %f\n",actualTime,energy);

  }while(actualTime < T);
                                                                                          
  return 0;
}


/*
 * La norma de un vector en dos dimensiones
 */
FLOAT norm( FLOAT vec1, FLOAT vec2 ){
  return sqrt(vec1*vec1 + vec2*vec2);
}
/*
 * El producto punto de dos vectores
 */
FLOAT dotP( FLOAT vec11, FLOAT vec12, FLOAT vec21, FLOAT vec22){
  return vec11*vec21 + vec21*vec22;
}
/*
 * El producto cruz de un vector 2D con el vector k
 */
FLOAT * crossP( FLOAT vec1, FLOAT vec2){
  FLOAT *vectemp = malloc(2*sizeof(FLOAT));
  FLOAT norma = norm(vec1,vec2);
  vectemp[0] = vec2/norma;
  vectemp[1] = -vec1/norma;
  return vectemp;
}
/*
 * Guarda memoria
 */
FLOAT * allocate( int n ){
  FLOAT * s;
  s = malloc( n*sizeof( FLOAT ) );
  return s;
}
/*
 * Inicializa las posiciones
 */
void initPos( FLOAT *xX, FLOAT *xY ){
  
  xX[0] = R*cos(0);
  xX[1] = R*cos(theta);
  xX[2] = R*cos(2*theta);
  xY[0] = R*sin(0);
  xY[1] = R*sin(theta);
  xY[2] = R*sin(2*theta);

}
/*
 * Inicializa las velocidades
 */
void initVel( FLOAT *vX, FLOAT *vY ){
 
  FLOAT x1 = R*cos(0);
  FLOAT x2 = R*cos(theta);
  FLOAT x3 = R*cos(2*theta);
  FLOAT y1 = R*sin(0);
  FLOAT y2 = R*sin(theta);
  FLOAT y3 = R*sin(2*theta);
  FLOAT *vt = crossP(x1,y1);
  vX[0] = v*vt[0];
  vY[0] = v*vt[1];
  vt = crossP(x2,y2);
  vX[1] = v*vt[0];
  vY[1] = v*vt[1];
  vt = crossP(x3,y3);
  vX[2] = v*vt[0];
  vY[2] = v*vt[1];

}
/*
 * Inicializa las masas
 */
void initMass( FLOAT *m ){
  m[0] = 1;
  m[1] = 2;
  m[2] = 3;
}

/*
 * Actualiza las aceleraciones
 */
void updateAcc( FLOAT *xX, FLOAT *xY, FLOAT *aX, FLOAT *aY , FLOAT *m){
  int i;
  int j;
  for( i = 0; i < n_masses; i++ ){
    aX[i] = 0;
    aY[i] = 0;
    for( j = 0; j < n_masses; j++ ){
      if( i != j ){
	FLOAT normi = pow(norm((xX[i] - xX[j]),(xY[i] - xY[j])), 3);
	FLOAT grav = (G*m[j])/normi;
	aX[i] += grav*(xX[j] - xX[i]);
	aY[i] += grav*(xY[j] - xY[i]);
      }
    }
  }
}

/*
 * Actualiza las posiciones
 */
void rungeKutta( FLOAT *xX, FLOAT *xY, FLOAT *vX, FLOAT *vY, FLOAT *m ,FLOAT *aX, FLOAT *aY, FLOAT dt){
  int i;
  // Vectores de pendientes intermedias Runge-Kutta
  FLOAT *k1x;
  FLOAT *k2x;
  FLOAT *k3x;
  FLOAT *k4x;
  FLOAT *kx_mean;
  FLOAT *k1y;
  FLOAT *k2y;
  FLOAT *k3y;
  FLOAT *k4y;
  FLOAT *ky_mean;
  FLOAT *k1vx;
  FLOAT *k2vx;
  FLOAT *k3vx;
  FLOAT *k4vx;
  FLOAT *kvx_mean;
  FLOAT *k1vy;
  FLOAT *k2vy; 
  FLOAT *k3vy;
  FLOAT *k4vy;
  FLOAT *kvy_mean;

  FLOAT *newXx;
  FLOAT *newXy;
  FLOAT *newVx;
  FLOAT *newVy;
  FLOAT *acx;
  FLOAT *acy;
  acx = allocate( n_masses );
  acy = allocate( n_masses );
  newXx = allocate( n_masses );
  newXy = allocate( n_masses );
  newVx = allocate( n_masses );
  newVy = allocate( n_masses );
 
  k1x = allocate( n_masses );
  k2x = allocate( n_masses );
  k3x = allocate( n_masses );
  k4x = allocate( n_masses );
  kx_mean = allocate( n_masses );
  k1y = allocate( n_masses );
  k2y = allocate( n_masses );
  k3y = allocate( n_masses );
  k4y = allocate( n_masses );
  ky_mean = allocate( n_masses );
  k1vx = allocate( n_masses );
  k2vx = allocate( n_masses );
  k3vx = allocate( n_masses );
  k4vx = allocate( n_masses );
  kvx_mean = allocate( n_masses );
  k1vy = allocate( n_masses );
  k2vy = allocate( n_masses );
  k3vy = allocate( n_masses );
  k4vy = allocate( n_masses ); 
  kvy_mean = allocate( n_masses ); 
  
  // UpdateAcc(xX,xY,aX,aY,m);
  for( i = 0; i < n_masses; i++){
    // Primera aproximación
    k1x[i] = vX[i];
    k1y[i] = vY[i];
    k1vx[i] = aX[i];
    k1vy[i] = aY[i];
  }
  
  for( i = 0; i < n_masses; i++){ 
    newXx[i] = xX[i] + k1x[i]*dt/2 ;
    newXy[i] = xY[i] + k1y[i]*dt/2 ;
    newVx[i] = vX[i] + k1vx[i]*dt/2;
    newVy[i] = vY[i] + k1vy[i]*dt/2;
    k2x[i] = newVx[i];
    k2y[i] = newVy[i];
  }
  
  updateAcc( newXx, newXy, acx, acy, m); 
  for( i = 0; i < n_masses; i++){ 
    k2vx[i] = acx[i];
    k2vy[i] = acy[i];
  }
  
  for( i = 0; i < n_masses; i++){ 
    // Tercera aproximación
    newXx[i] = newXx[i] + k2x[i]*dt/2;
    newXy[i] = newXy[i] + k2y[i]*dt/2;
    newVx[i] = newVx[i] + k2vx[i]*dt/2;
    newVy[i] = newVy[i] + k2vy[i]*dt/2;
    k3x[i] = newVx[i] ;
    k3y[i] = newVy[i] ;
  }
  
  updateAcc( newXx, newXy, acx, acy, m); 
  for( i = 0; i < n_masses; i++){ 
    k3vx[i] = acx[i];
    k3vy[i] = acy[i];
  }

  for( i = 0; i < n_masses; i++){ 
    // Cuarta aproximación
    newXx[i] = newXx[i] + k3x[i]*dt;
    newXy[i] = newXy[i] + k3y[i]*dt;
    newVx[i] = newVx[i] + k3vx[i]*dt;
    newVy[i] = newVy[i] + k3vy[i]*dt;
    k4x[i] = newVx[i];
    k4y[i] = newVy[i];
  }
  updateAcc( newXx, newXy, acx, acy, m); 
  for( i = 0; i < n_masses; i++){ 
    k4vx[i] = acx[i];
    k4vy[i] = acy[i];
  }

  for( i = 0; i < n_masses; i++){ 
    // Promedio
    kx_mean[i] = (1/6.0)*(k1x[i] + 2*k2x[i] + 2*k3x[i] + k4x[i]);
    ky_mean[i] = (1/6.0)*(k1y[i] + 2*k2y[i] + 2*k3y[i] + k4y[i]);
    kvx_mean[i] = (1/6.0)*(k1vx[i] + 2*k2vx[i] + 2*k3vx[i] + k4vx[i]);
    kvy_mean[i] = (1/6.0)*(k1vy[i] + 2*k2vy[i] + 2*k3vy[i] + k4vy[i]);
    // Finalmente, los valores actualizados
    xX[i] = xX[i] + kx_mean[i]*dt;
    xY[i] = xY[i] + ky_mean[i]*dt;
    vX[i] = vX[i] + kvx_mean[i]*dt;
    vY[i] = vY[i] + kvy_mean[i]*dt;
  }
}
/*
 * Calcula la energía total del sistema
 */
FLOAT calcEnergy(FLOAT *xX, FLOAT *xY, FLOAT *vx, FLOAT *vy, FLOAT *m){
  FLOAT ener = 0;
  int i,j;
  for( i = 0; i < n_masses; i++){
    ener += (1/2.0)*vx[i]*vx[i];
    ener += (1/2.0)*vy[i]*vy[i];
    for( j = 0; j < n_masses; j++){
      if( i != j){
	FLOAT normi = norm((xX[i] - xX[j]),(xY[i] - xY[j]));
	ener += -G*m[i]*m[j]/normi;
      }
    } 
  }
  return ener;
}



