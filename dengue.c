/* Codigo para estudiar enfermedades vectoriales (particularmente dengue)

Names
poblations: *) V: mosquitoes (V=W+Z)
            *) W: susceptibles mosquitoes
            *) Z: infected mosquitoes
            *) N: humans
            *) S: susceptibles humans
            *) I: infected humans
            *) R: recovered humans

events: *) lamda=mu: birth and deth rate of moquitoes
        *) bitting_rate: bitting rate of mosquitoes
        *) Trasn_HV: probability to transfer the virus form humans to mosquitoes (given that a susceptible mosquito has bitten to an infected human)
        *) Trasn_VH: probability to transfer the virus form mosquitoes to humans (given that a infected mosquito has bitten to an susceptible human)
        *) IncubationPeriod: incubation period of the virus inside the human. Time while the human can transmit the virus to vectors (before this time is not possible the virus transmission)
        *) lamdaH = mH : death and birth rate of humans

Associated equations
dW/dt=lamda*V-bitting_rate*Trans_HV*I*W/N-lamda*W
dZ/dt= bitting_rate*Trans_HV*I*W/N-lamda*Z

dS/dt= -bitting_rate*Trans_VH*Z*S/N + lamdaH*N - lamdaH*S
dI/dt= bitting_rate*Trans_VH*Z*S/N - (1/IncubationPeriod)*I - lamdaH*I
dR/dt= (1/IncubationPeriod)*I - lamdaH*R
*/ 

#include "dengue.h"
#include "algoritmos.c" /* generation of random numbers */

/* The global variables are defined*/
common C; /* I should put it as a local varible structure*/

long pobla[LIMfilas][LIMcolumnas][POBLACIONES];/* population array: two first dimension are the spatial index */
long neventos[LIMfilas][LIMcolumnas][EVENTOS+14]; /*population array: two first dimension are the spatial index, the + is asociated to the vuelos to different units*/
double Lbd[LIMfilas][LIMcolumnas][EVENTOS];

#include "rates_dengue.c" /* rates that are not constant in te time */
#include "f-auxiliares_dengue.c" /* calculo de coeficientes y rutinas de entrada */

#include "vm_dengue.c" /* rutina de calculo de valores medios de poblaciones */
#include "rk2.c" /* Runge-Kutta integrator */
#include "deriv_dengue.c"
#include "reparto_dengue.c"


/* Cuerpo principal del programa */

int main()
{
   long idum= 3480;  /* Seed to generate pseudorandom numbers*/
   
   double t=0.0, *pointerLbd;

/* POPULATIONS
   0 susceptibles mosquitoes
   1 infected mosquitoes
   2 susceptibles humans
   3 infected humans
   4 recovered humans */

/* EVENTS
   0 -- susceptibles mosquitoes mortality
   1 -- infection of susceptibles mosquitoes
   2 -- infected mosquitoes mortality
   3 -- mosquitos birth
   4 -- susceptibles humans mortality
   5 -- infection of susceptibles humans
   6 -- infected humans mortality
   7 -- recovery of infected humans
   8 -- recovered humans mortality
   9 -- humans birth
   10 - flight of susceptibles mosquitoes
   11 - flight of infected mosquitoes
*/

   char *poblaname[POBLACIONES], mastername[40], *eventoname[EVENTOS+14]; 
   char poblaN[POBLACIONES][48], eventoN[EVENTOS+14][48]; 
   FILE *FilePobla[POBLACIONES], *FileEventos[EVENTOS+14]; 
/*lectura de datos e impresion de resultados*/

   int i,l,o,Toffset; 
   long *pointereventos;
/* tiempo del transitorio y de las simulacion (en días)  a ser leidos. fecha de inicio queda fija en el 1 de julio */

   int Transitorio=0, Tsimulado=0, TransitorioLeido=-1, repite;
   float Dt=1./PASO; /* paso de tiempo en fraccion dias */

   int Plist[POBLACIONES], Elist[EVENTOS+14]; 
/* listado de eventos y poblaciones para imprimir -Print- */

   float temperaturavec[5000]; /* vector with the temperature values*/
   float Temp=0.0;
   int tiempo_inicial=0; /* tiempo en dias al comenzar el cálculo */

   pointerLbd =Lbd[0][0]; /* pointer para recorrer los lambdas */
   pointereventos=neventos[0][0]; /* pointer para recorrer eventos */

/* asigno los pointers de los nombres POCO ELEGANTE */
   for (i=0;i<POBLACIONES; i++) 
      poblaname[i]=poblaN[i];
   for (i=0;i<(EVENTOS+14); i++) 
      eventoname[i]=eventoN[i];

/*lee los datos de entrada y abre los files, todos los argumentos son salida */

   if( (i=getdata(FilePobla, FileEventos, Plist, Elist, mastername, poblaname, eventoname, &Transitorio, &Tsimulado, &idum, temperaturavec)) < 0 ) {
      fprintf(stderr,"Error al leer %s es %d\n", DENGUEIMPUT,i);};

   fprintf(stderr,"Transitorio %d  Simulado %d semilla %ld \n", Transitorio, Tsimulado, idum);

   TransitorioLeido=Transitorio;

/* inicializo las constantes de evolucion que no dependen de la Temperatura */
#include "CONSTANTES_dengue"

/* REPETICION para la estadistica */
   for (repite=1; repite <= REPITE; repite++) {
      fprintf(stderr, "Repetición %d\n",repite); 
/* inicializo las poblaciones */
     pobla_iniciales(&Transitorio, &tiempo_inicial);//,&idum); 
    

     for(i=0; i< (EVENTOS+14) * LIMfilas * LIMcolumnas; i++) 
	*(pointereventos+i)=0; /* pongo en cero los eventos del dia */
/* inicio el transitorio */
      for (l=tiempo_inicial; l< Transitorio+tiempo_inicial ; l++) {  
         Temp=temperaturavec[l-TsaveCI]; /* doy valores a Temp */
         coeficientes(Temp); /* llamo a la funcion coeficientes que me da los valores para los rates que dependen de la temperatura*/

	 if((l == 0) && (tiempo_inicial != 0)) save_transitorio(pobla);
	 /* corre un día */
/* Hernan pone a cero los eventos */
     for(i=0; i< (EVENTOS+14) * LIMfilas * LIMcolumnas; i++) 
	*(pointereventos+i)=0; /* pongo en cero los eventos del dia */
         for(o=0; o< PASO; o++){
            t=l+(o+1.)*Dt;
	    for(i=0; i< TOTALEVENTOS; i++)
	       *(pointerLbd+i)=0; /* pongo a cero el lambda */
	    rk2(deri, Lbd[0][0], TOTALEVENTOS, t, Dt);
	    updatepobla(&idum,t);  
	 };
      };

      if(Transitorio+tiempo_inicial-TsaveCI > 0) {
         Temp=temperaturavec[Transitorio+tiempo_inicial-TsaveCI-1];  /*miro la temperatura en el vector al final del transitorio*/
      }
      else {
         Temp=temperaturavec[0];
      }     

      for (i=0; (i<POBLACIONES) && (Plist[i] >=0 ); i++)
         fprintf(FilePobla[i],"%d\n", repite); 
      for (i=0; (i<EVENTOS) && (Elist[i] >=0 ); i++)
  	 fprintf(FileEventos[i],"%d\n",repite);
/* Pongo un cartel en los files */
/* SIMULACION */
      if((l == 0) && (tiempo_inicial != 0) ) save_transitorio(pobla);  /*guardo las cond iniciales*/

/* Hernan corrige la salida de Y */
      Toffset= l - (TransitorioLeido+TsaveCI);
      savedata(FilePobla,FileEventos,Plist,Elist, t=(float)Toffset); /* guardo la salida */ 

/* loop de tiempo despues de guardar la CI */
      for (l=tiempo_inicial+Transitorio; l< -1+Tsimulado+Transitorio+tiempo_inicial ; l++) { 
         Temp=temperaturavec[l-TsaveCI]; /* doy valores a Temp */
         coeficientes(Temp); /* llamo a la funcion coeficientes que me da los valores para los rates que dependen de la temperatura*/
      
         for(i=0; i< (EVENTOS+14) * LIMfilas * LIMcolumnas; i++) 
	    *(pointereventos+i)=0; /* pongo en cero los eventos del dia */
	 for(o=0; o< PASO; o++) {
 	    t=l+o*Dt;
            
	    for(i=0; i< TOTALEVENTOS; i++)
	       *(pointerLbd+i)=0; /* pongo a cero el lambda */
         
	    rk2(deri, Lbd[0][0], TOTALEVENTOS, t, Dt);
	    updatepobla(&idum,t);



         }; /*fin del for del indice o*/
      
        Toffset= l - (TransitorioLeido+TsaveCI) + 1;


      savedata(FilePobla,FileEventos,Plist,Elist, t=(float)Toffset); /*guardo las pob, eventos y temp */
            
      };  /*fin del loop del tiempo, indice l*/      
      
     
      Transitorio= TransitorioLeido; /* Transitorio va a ser decrementado en TsaveCI al leer el archivo de CI generado */
   
   

   }/* fin del loop de repeticion */

   nuevasemilla(-idum); /* guarde la semilla el menos para ran1 de nrc con + para ran33 */

   exit(0);
}
