#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#define PI 3.14159265358979323846

/*
 * La norma de un vector en dos dimensiones
 */
double norm( double vec_1, double vec_2 );
/*
 * El producto punto de dos vectores
 */
double dotP( double vec_1_1, double vec_1_2, double vec_2_1, double vec_2_2);
/*
 * El producto cruz de un vector 2D con el vector k
 */
double* crossP( double* vec);

int main(int argc, char **argv){
 
  //------------------------------------------------------------------------
  // Constantes
  //------------------------------------------------------------------------

  // En unidades de masa solar, 1au, años
  double solar_mass = 1.0;
  double m1 = 1.0;
  double m2 = 2.0;
  double m3 = 3.0;
  double G = 39.486;
  double dt =  (double) atof(argv[1]);
  double T = 5000.0;

  // Archivos
  FILE *datax1;
  FILE *datax2;
  FILE *datax3;
  FILE *datav1;
  FILE *datav2;
  FILE *datav3;
  FILE *energyData;
  datax1  = fopen("x1.data", "w");
  datax2 = fopen("x2.data", "w");
  datax3 = fopen("x3.data", "w");
  datav1 = fopen("v1.data", "w");
  datav2 = fopen("v2.data", "w");
  datav3 = fopen("v3.data", "w");
  energyData = fopen("energy.data", "w");

  // Condiciones iniciales
  double R = 100.0;
  double r = sqrt(3.0)*R;
  double theta = PI/3.0;
  double v = sqrt((11.0*G*solar_mass)/(3.0*r));
  double n_dim = 2;
  
  //------------------------------------------------------------------------
  // Variables 
  //------------------------------------------------------------------------
    
  double* x1; 
  double* x2;
  double* x3;
  double* v1;
  double* v2;
  double* v3;
  x1 = malloc(n_dim*sizeof(double));
  x2 = malloc(n_dim*sizeof(double));
  x3 = malloc(n_dim*sizeof(double));
  v1 = malloc(n_dim*sizeof(double));
  v2 = malloc(n_dim*sizeof(double));
  v3 = malloc(n_dim*sizeof(double));
    
  //------------------------------------------------------------------------
  // Valores iniciales
  //------------------------------------------------------------------------
  
  x1[0] = 0;
  x1[1] = R;
  x2[0] = R*cos(theta);
  x2[1] = -R*sin(theta);
  x3[0] = -R*cos(theta);
  x3[1] = -R*sin(theta);

  double* tempVec;
  tempVec = crossP( x1 );
  v1[0] = (v/R)*(*tempVec[0]);
  v1[1] = (v/R)*(*tempVec[1]);
  tempVec = crossP( x2 );
  v2[0] = (v/R)*(*tempVec[0]);
  v2[1] = (v/R)*(*tempVec[1]);
  tempVec = crossP( x3 );
  v3[0] = (v/R)*(*tempVec[0]);
  v3[1] = (v/R)*(*tempVec[1]);

  //Aqui va un while

  //------------------------------------------------------------------------
  // 1st second order DE Runge - Kutta (x1,v1)
  //------------------------------------------------------------------------
  
  double* tempx1;
  double* tempv1;
  double tempa1;
  double tempa2;  
  double norm1 = norm(x1[0] - x2[0], x1[1] - x2[1]);
  double norm2 = norm(x1[0] - x3[0], x1[1] - x3[1]);
  tempa1 = ((G*m2/(norm1*norm1*norm1))*(x2[0]-x1[0])) +  ((G*m3/(norm2*norm2*norm2))*(x3[0]-x1[0]));
  tempa2 = ((G*m2/(norm1*norm1*norm1))*(x2[1]-x1[1])) +  ((G*m3/(norm2*norm2*norm2))*(x3[1]-x1[1]));

  // runge-kutta step

  
  
                                                                                          
  return 0;
}


double norm( double vec_1, double vec_2 ){
  return sqrt(vec_1*vec_1 + vec_2*vec_2);
}

double dotP( double vec_1_1, double vec_1_2, double vec_2_1, double vec_2_2){
  return vec_1_1*vec_2_1 + vec_1_2*vec_2_2;
}

double* crossP( double* vec){
  double* newVec;
  newVec = malloc(2*sizeof(double));
  newVec[0] = vec[1];
  newVec[1] = -vec[0];
  return newVec;
}
