/* Calculating the rates derivates */
void deri(int dimension, double l[], double dl[], double t)
{/* l is the lambda in a vector way, dl is the derivate in a vector way and t the time that here we don't used. */
double Tnl[LIMfilas][LIMcolumnas][EVENTOS],vm[LIMfilas][LIMcolumnas][POBLACIONES];
int i,j,m;

float MosquitosTotales, VectorMortality;
long H0;


VM(vm,Tnl); /* function with mean value of the population for each time step */
            /* Tnl is the non-linear term --> Taylor first order factor*/

for (i=0; i < LIMfilas; i++)
	for (j=0; j < LIMcolumnas; j++)
	   {    
                
                VectorMortality=C.mm/C.AmbientCapacity[i][j];
                MosquitosTotales=vm[i][j][0]+vm[i][j][1]; 
                H0 = pobla[i][j][2]+pobla[i][j][3]+pobla[i][j][4];      

                m=(LIMcolumnas *i+ j)*EVENTOS;
		(dl+m)[0] = VectorMortality * (MosquitosTotales * vm[i][j][0] + Tnl[i][j][0]);                                                 /* susceptibles mosquitoes mortality */
		(dl+m)[1] = (rateInfection_ordenCero(H0,pobla[i][j][3]) + rateInfection_PrimerOrden(H0,Tnl[i][j][1]))*vm[i][j][0]; 	       /* infection of susceptibles mosquitoes */
		(dl+m)[2] = VectorMortality * (MosquitosTotales * vm[i][j][1] + Tnl[i][j][2]);                                                 /* infected mosquitoes mortality */
		(dl+m)[3] = C.nm * (vm[i][j][0]+vm[i][j][1]); 	                                                                               /* mosquitoes birth */
		(dl+m)[4] = C.mh * vm[i][j][2];                                                                                                /* susceptibles humans mortality */
		(dl+m)[5] = (rateInfection_ordenCero(H0,pobla[i][j][1]) + rateInfection_PrimerOrden(H0,Tnl[i][j][5]))*vm[i][j][2]; 	       /* infection of susceptibles humans */
                (dl+m)[6] = C.mh * vm[i][j][3];                                                                                                /* infected humans mortality*/
                (dl+m)[7] = (1/C.IncubationPeriod)*vm[i][j][3];                                                                                /* recovery of infected humans */
		(dl+m)[8] = C.mh * vm[i][j][4];                                                                                                /* recovered humans mortality */
		(dl+m)[9] = C.nh*(vm[i][j][2]+vm[i][j][3]+vm[i][j][4]);                                                                        /* humans birth */
                (dl+m)[10] = C.vuelo[i][j]*vm[i][j][0];                                                                                        /* flight of susceptibles mosquitoes */
                (dl+m)[11] = C.vuelo[i][j]*vm[i][j][1];                                                                                        /* flight of infected mosquitoes */
         
	         
         
		};
}



