#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define DENGUEIMPUT "dengue.inp" /* imput file for the code */
#define LIMfilas 21 /*total number of rows in the grid  */ 
#define LIMcolumnas 21 /*total number of columns in the grid  */
#define PASO 12
#define REPITE 100 /* Number of repetitions */
#define TsaveCI -0 /* Transition time (time before the simulation) */

typedef struct Common
	{float mm,nm,bitr,trans_mh,trans_hm,IncubationPeriod;
         float AmbientCapacity[LIMfilas][LIMcolumnas],gente[LIMfilas][LIMcolumnas],vuelo[LIMfilas][LIMcolumnas];
         double pvuel[8];
         float mh, nh;
         float V0,lambda; 
         int centro,Inf0, Nh, TotalVectors;
  	} common;  /* define el type de common  for the Common structure. This contain pointers to the variables 
that are used in subrutines such as the calculation of the derivarives*/ /* Some of this variables are in CONSTANTES*/


#define JULIO 184         /* off set del 1 de julio a 1 de enero */ 
#define EVENTOS 12        /*number of events per unit*/
#define POBLACIONES  5	  /* number of populations per unit */
#define TOTALPOBLA LIMcolumnas * LIMfilas * POBLACIONES /* total of populations in the grid */
#define TOTALEVENTOS   LIMcolumnas * LIMfilas * EVENTOS /* total of events in the grid */

/* Prototipos de funciones */
float ran1(long *);         /* generates a random number with a uniform distribution*/
float gammln(float );       /* Ln gamma distribution. It is necessary for the Poisson distribution*/
int poidev(float , long *); /* Poisson distribution*/
float gasdev(long *);       /* generates random numbers with a normal N(0,1) distribution*/

/* simple multinomial */
unsigned long BINV(double, long ,long *);
unsigned long BTRD(double ,long ,long *);
unsigned long binomial(double ,long ,  long *);
void multinomialS(int , int , int ,double [],int [],long *);
void multinomialR(int ,int ,int ,double [],int [],long *);



