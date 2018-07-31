/* We make the allocation of population increases */

void reparto (int i, int j, double lam[],long  paux[][LIMcolumnas][POBLACIONES], int count[], long *idum, double t)
{int N1=0, N2=0, N3=0;
int a,e,im, jm, ip, jp;
//int cnt[EVENTOS]={0};
int cnt[EVENTOS];
double prob;
int count_aux[2];

/*the allocation of the positive events */
			/* mosquitos suceptibles */
                        count_aux[0]=0;
			if(lam[0] > 0) /* careful with the division by zero */
				{prob=Lbd[i][j][0]/lam[0];

				N1=binomial(prob, count[0], idum); /*death*/
				N2=count[0]-N1; /*they become infected or fly to another adjacent cell*/
                                neventos[i][j][0]+= N1; /* mortality susceptibles mosquitoes events*/
                                
                                if(N2>0) {
                                   prob=Lbd[i][j][1]/N2;
                                   count_aux[0]=N2-binomial(prob,N2,idum);/* they fly */
                                   N2-=count_aux[0]; /* they become infected*/
                                }
				paux[i][j][1] += N2 ; /*adding to  the infected mosquitoes*/
				neventos[i][j][1]+= N2; /* passage events from susceptibles mosquitoes to infected*/
				};

                        /*infected mosquitoes*/
                        count_aux[1]=0;
			if(lam[1] > 0)
				{prob=Lbd[i][j][2]/lam[1];
				N1=binomial(prob, count[1], idum); /* death*/
                                count_aux[1]=count[1]-N1;
				neventos[i][j][2]+= N1;
                                
				/*if(N2 != 0)
					{neventos[i][j][3]+= N2;
					fprintf(stderr,"This should be zero N2=%d i=%d j=%d\n",N2,i,j);
					}*/
				}
                        /*susceptibles humans*/
			if(lam[2] > 0)
				{prob=Lbd[i][j][4]/lam[2];
				N1=binomial(prob, count[2], idum); /* death*/
				N2=count[2]-N1; /* infections */
				paux[i][j][3] += N2 ; /* adding to the infected */
				neventos[i][j][4]+= N1; /* adding to the death event */
				neventos[i][j][5]+= N2; /* adding to the infection event*/
				}
			/* infected humans */
			if(lam[3] > 0)
				{prob=Lbd[i][j][6]/lam[3];
				N1=binomial(prob, count[3], idum); /* death*/
				N2=count[3]-N1; /* recoveries */
				paux[i][j][4] += N2 ; /* adding to the recovered population */
				neventos[i][j][6]+= N1; /* adding to the death event */
				neventos[i][j][7]+= N2; /* adding to the recovery event */
				}
			/* recovered humans */
			if(lam[4] > 0)
				{//prob=Lbd[i][j][8]/lam[4];
				//N1=binomial(prob, count[4], idum); // mueren
                                N1=count[4]; /*death*/
				//N2=count[4]-N1; /* this should be zero because this population has one "negative" event*/
				neventos[i][j][8]+= N1; /*adding to the death event*/
				}
                        /* the births are added */
                        /* for mosquitos*/
                        N1=count[5];
                        paux[i][j][0]+=N1;
			neventos[i][j][3]+=N1;
                        /* for humans */
                        N1=count[6];
                        paux[i][j][2]+=N1;
                        neventos[i][j][9]+=N1;

                       
                        /* finally the flight */
                        im=i-1; jm=j-1; ip=i+1;jp=j+1;
       			if( (im <= 0) ) im=0;
                	if( (jm <= 0) ) jm=0;
                	if( (jp >= LIMcolumnas) ) jp=LIMcolumnas-1;
                        if( (ip >= LIMfilas) ) ip=LIMfilas-1;
                        for(a=0;a<2;a++){
                           for(e=0;e<EVENTOS;e++)
                             cnt[e]=0;
                           if(count_aux[a]>0){
                              multinomial(0,7,count_aux[a],C.pvuel,cnt,idum);
                              paux[im][jm][a]+= cnt[0]; 
			      paux[im][j][a]+= cnt[1];
			      paux[im][jp][a]+= cnt[2];
			      paux[i][jp][a]+= cnt[3];
			      paux[ip][jp][a]+= cnt[4];
			      paux[ip][j][a]+= cnt[5];
			      paux[ip][jm][a]+= cnt[6];
			      paux[i][jm][a]+= cnt[7];
                              for(e=0;e<8;e++){
                                 neventos[i][j][(10+e+8*a)]+=cnt[e];
                                 if(a==1)
                                    neventos[i][j][25]+=cnt[e];
                                 if(a==0)
                                    neventos[i][j][10]+=cnt[e];
                              }
                           }
                        }                
}



void updatepobla (long *idum,float t)
{
  int i,j,k,e, count[EVENTOS]={0};
  long Etotales=0;
  long resto[LIMfilas][LIMcolumnas][POBLACIONES], restototal=0, oldresto=-1;
  long paux[LIMfilas][LIMcolumnas][POBLACIONES];
  double Ltotal, MAXl=0.0, lam[EVENTOS];


/* A copy of the population is made */ 
  for (i=0; i < LIMfilas; i++)
    for (j=0; j< LIMcolumnas; j++){
      
      for(k=0; k<POBLACIONES ; k++)
        paux[i][j][k]=pobla[i][j][k];
      
    }



  for (i=0; i< LIMfilas; i++)
    for (j=0; j< LIMcolumnas; j++) {
      /* the count variable is set to zero */
      for(e=0;e<EVENTOS;e++)
         count[e]=0;

      Ltotal=0.0;
      for (k=0; k< (POBLACIONES); k++) {
        lam[k] = Lbd[i][j][2*k]+Lbd[i][j][2*k+1];
	resto[i][j][k]=0;
	Ltotal += lam[k];
/* lam[k] has the sum of the lambdas that subtract to the population k. As it is indicated in the main*/

	MAXl= (MAXl >= Ltotal? MAXl : Ltotal); 
      }
      
      lam[0] +=Lbd[i][j][10]; /* flights */
      lam[1] +=(Lbd[i][j][11]-Lbd[i][j][3]);
      lam[4] -= Lbd[i][j][9]; /* births */
      Ltotal +=(Lbd[i][j][10]+Lbd[i][j][11]-Lbd[i][j][3]-Lbd[i][j][9]);
      
      Etotales= poidev(Ltotal, idum); /* total number the events in a unit */
      multinomial(0, POBLACIONES-1, Etotales, lam, count, idum);

      count[5]=poidev(Lbd[i][j][3], idum); /*has the negative events for a non-existen population, but it is positive por the susceptibles mosquitoes */
      count[6]=poidev(Lbd[i][j][9], idum); /* has the negative events for a non-existen population, but it is positive por the susceptibles humans */
      
     
/* in count variable we have the negative events for each population. The values of count come from the multinomial */
      for (k=0; k< POBLACIONES; k++) {
	if (count[k] > pobla[i][j][k]) {/* temporary excess of events */
	  resto[i][j][k]= count[k]-pobla[i][j][k];
	  count[k]= pobla[i][j][k];
	  paux[i][j][k]=0; /* the population is set to zero */
	  restototal += resto[i][j][k]; 
	}
	else /* subtraction */ {
          if(paux[i][j][k]<0)
             printf("WARNING\n");
          paux[i][j][k] -= count[k];
	};
      };

      reparto(i, j, lam, paux, count, idum,t);
    }; //end of j and i loops

/* population update */
  for (i=0; i< LIMfilas; i++)
    for (j=0; j< LIMcolumnas; j++){
      for (k=0; k< POBLACIONES; k++)
	pobla[i][j][k]=paux[i][j][k];
}
  while ((restototal != oldresto) && (restototal > 0)) {
    oldresto=restototal; restototal=0;
    for (i=0; i< LIMfilas; i++)
      for (j=0; j< LIMcolumnas; j++) {
        for(k=0; k < (POBLACIONES) ; k++) {
          lam[k] = Lbd[i][j][2*k]+Lbd[i][j][2*k+1];
          count[k]=resto[i][j][k];
        };
        lam[0] +=Lbd[i][j][10];
        lam[1] +=(Lbd[i][j][11]-Lbd[i][j][3]);
        lam[4] -= Lbd[i][j][9]; /* births */
        for(k=0; k < POBLACIONES ; k++) {
          paux[i][j][k]=pobla[i][j][k];
          if (count[k] > pobla[i][j][k]) {
            resto[i][j][k]= count[k]-pobla[i][j][k];
	    count[k]= pobla[i][j][k];
	    paux[i][j][k]=0;
	    restototal += resto[i][j][k];
          }
          else {
            paux[i][j][k] -= count[k];
	  };
        }; 

        reparto(i, j, lam, paux, count, idum,t);

      };
      
      for (i=0; i< LIMfilas; i++)
        for (j=0; j< LIMcolumnas; j++)
	  for (k=0; k< POBLACIONES; k++)
	    pobla[i][j][k]=paux[i][j][k];
  }; /* while loop */

}/*End update*/



