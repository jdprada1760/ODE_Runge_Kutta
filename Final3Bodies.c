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
#define v sqrt((1.10*G*solar_mass)/(3.0*r))

FLOAT norm( FLOAT vec1, FLOAT vec2 );
FLOAT dotP( FLOAT vec11, FLOAT vec12, FLOAT vec21, FLOAT vec22);
FLOAT * crossP( FLOAT vec1, FLOAT vec2);
FLOAT * allocate( int n );
void initPos( FLOAT *xX, FLOAT *xY );
void initVel( FLOAT *vX, FLOAT *vY );
void initMass( FLOAT *m );
void updateAcc( FLOAT *xX, FLOAT *xY, FLOAT *aX, FLOAT *aY , FLOAT *m);
FLOAT getAccel( FLOAT *xX, FLOAT *xY, FLOAT *m, FLOAT newXx, FLOAT newYy, int j);
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
 * Obtiene la aceleración del elemento iésimo en función de las posiciones relativas
 */ 
FLOAT getAccel( FLOAT *xX, FLOAT *xY, FLOAT *m, FLOAT newXx, FLOAT newYy, int j){
  int i;
  FLOAT ans = 0;
  for( i = 0; (i < n_masses); i++){
    if( i != j ){
      FLOAT normi = pow(norm((xX[i] - newXx),(xY[i] - newYy)), 3);
      ans += (G*m[j]*(xX[i]-newXx))/normi;
    }
  }
  return ans;
}

/*
 * Actualiza las posiciones
 */
void rungeKutta( FLOAT *xX, FLOAT *xY, FLOAT *vX, FLOAT *vY, FLOAT *m ,FLOAT *aX, FLOAT *aY, FLOAT dt){
  int i;
  //updateAcc(xX,xY,aX,aY,m);
  for( i = 0; i < n_masses; i++){
    
    // Primera aproximación
    FLOAT k1x = vX[i];
    FLOAT k1y = vY[i];
    FLOAT k1vx = aX[i];
    FLOAT k1vy = aY[i];
    
    // Segunda aproximación
    FLOAT newXx = xX[i] + k1x*dt/2 ;
    FLOAT newXy = xY[i] + k1y*dt/2 ;
    FLOAT newVx = vX[i] + k1vx*dt/2;
    FLOAT newVy = vY[i] + k1vy*dt/2;
    FLOAT k2x = newVx ;
    FLOAT k2y = newVy ;
    FLOAT k2vx = getAccel(xX,xY,m,newXx,newXy,i);
    FLOAT k2vy = getAccel(xY,xX,m,newXy,newXx,i);
    
    // Tercera aproximación
    newXx = newXx + k2x*dt/2;
    newXy = newXy + k2y*dt/2;
    newXx = newVx + k2vx*dt/2;
    newXx = newVy + k2vy*dt/2;
    FLOAT k3x = newVx ;
    FLOAT k3y = newVy ;
    FLOAT k3vx = getAccel(xX,xY,m,newXx,newXy,i);
    FLOAT k3vy = getAccel(xY,xX,m,newXy,newXx,i);
    
    // Cuarta aproximación
    newXx = newXx + k3x*dt;
    newXy = newXy + k3y*dt;
    newXx = newVx + k3vx*dt;
    newXx = newVy + k3vy*dt;
    FLOAT k4x = newVx ;
    FLOAT k4y = newVy ;
    FLOAT k4vx = getAccel(xX,xY,m,newXx,newXy,i);
    FLOAT k4vy = getAccel(xY,xX,m,newXy,newXx,i);

    // Promedio
    FLOAT kx_mean = (1/6.0)*(k1x + 2*k2x + 2*k3x + k4x);
    FLOAT ky_mean = (1/6.0)*(k1y + 2*k2y + 2*k3y + k4y);
    FLOAT kvx_mean = (1/6.0)*(k1vx + 2*k2vx + 2*k3vx + k4vx);
    FLOAT kvy_mean = (1/6.0)*(k1vy + 2*k2vy + 2*k3vy + k4vy);

    // Finalmente, los valores actualizados
    xX[i] = xX[i] + kx_mean*dt;
    xY[i] = xY[i] + ky_mean*dt;
    vX[i] = vX[i] + kvx_mean*dt;
    vY[i] = vY[i] + kvy_mean*dt;

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
	ener += G*m[i]*m[j]/normi;
      }
    } 
  }
  return ener;
}



