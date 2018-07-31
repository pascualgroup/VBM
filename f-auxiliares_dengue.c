/* save the data after a day */
void savedata(FILE *FilePobla[], FILE *FileEventos[], int Plist[], int Elist[], float time)	
{
int err=0,i,j;

while ((Plist[err] >= 0) && (err < (POBLACIONES)) )
	{ /* type the populations */
 	fprintf(FilePobla[err], "%f \n", time);
	for (i=0; i< LIMfilas ; i++)
		{for (j=0; j< LIMcolumnas; j++)
			fprintf(FilePobla[err], "%ld\t",pobla[i][j][err]);
		fprintf(FilePobla[err], "\n"); 
		}
	err++ ; /* advances the index */
	}


err=0;  /* type the events */
while ((Elist[err] >= 0) && (err < EVENTOS+14) ) 
	{
	fprintf(FileEventos[err], "%f \n", time);
	for (i=0; i< LIMfilas ; i++)
		{for (j=0; j< LIMcolumnas; j++)
			fprintf(FileEventos[err], "%ld\t",neventos[i][j][err]);
		fprintf(FileEventos[err], "\n");
		}
	err++ ; 
	} /* advances the index */
}

int getdata(FILE *FilePobla[], FILE *FileEventos[], int Plist[], int Elist[],
	char *mastername, char *poblaname[], char *eventoname[],//
	int *Transitorio, int *Tsimulado, long *idum, float temperaturavec[])
{
int Npobla=0, Neventos=0;
int err=0,i;
char aux2[10]="NADA", aux1[50]="NADA";
FILE *inputfile;
inputfile=fopen(DENGUEIMPUT,"r");
if(inputfile == (FILE *) NULL)
	return -1; /* error of reading the file */
fscanf(inputfile,"%s\n%d %d",mastername,&Npobla,&Neventos);
fprintf(stderr,"Lei %s\n%d %d\n",mastername,Npobla,Neventos);

/* Open the population files */
if ((Npobla < 0) || (Npobla > POBLACIONES) )
	return err=-2;
else
	{for (i=0; i< Npobla; i++)
		fscanf(inputfile,"%d", &Plist[i]);
	} 
	for (i=Npobla; i< POBLACIONES; i++)
		Plist[i]=-1;

for (i=0; i < Npobla ; i++) 
	{
	strcpy(aux1,mastername);/* copy the master name to the auxiliary */
	
	sprintf(aux2,"-P%d.dat",Plist[i]); /* the indicator is added */
	strcat(aux1,aux2); /*the indicator is added */
	strcpy(poblaname[i],aux1); /* save the name */
	FilePobla[i]=fopen(poblaname[i],"w"); /* pass the pointer to the filepointers vector */
	};


/* the event section */
if ((Neventos < 0) || (Neventos > (EVENTOS+14)) ) 
	return err=-3;
else
	{for (i=0; i< Neventos; i++)
		fscanf(inputfile,"%d", (Elist+i));
	}
	for (i=Neventos; i< (EVENTOS+14); i++) 
		Elist[i]=-1;


for (i=0; i < Neventos ; i++) 
	{
	strcpy(aux1,mastername);/* copy the master name to the auxiliary */
	
	sprintf(aux2,"-E%d.dat",Elist[i]);
	strcat(aux1, aux2); /* the indicator is added */
	strcpy(eventoname[i],aux1); /* save the name */
	FileEventos[i]=fopen(eventoname[i],"w"); /* pass the pointer to the filepointers vector*/
	};
fscanf(inputfile,"%d %d %i", Transitorio, Tsimulado, &C.Nh); /* C.Nh is the total number of humans in the grid */
fclose(inputfile);

if ((inputfile=fopen("semilla.inp","r") )!= NULL)
	fscanf(inputfile, "%ld", idum);
else
	{
	printf("preciso semilla ");
	scanf("%ld", idum);
	}

int filatemp;
FILE *tempylluv;  /*The file in one column*/
tempylluv = fopen("Temp.dat", "r");
if(tempylluv == (FILE *) NULL)
   return -4; /* error of reading the file */
for (filatemp=0; filatemp < *Tsimulado + *Transitorio; filatemp++) {
   if(EOF == fscanf(tempylluv,"%f\n",&temperaturavec[filatemp])) {fprintf(stderr,"NO MORE WEATHER DATA\n"); exit(-1);}
}
fclose(tempylluv);
return err;

}


void save_transitorio(long p[][LIMcolumnas][POBLACIONES]) {
/* p is the pointer to de poblational vector */
   
   FILE *tedio;
   int i,j,k;

   if((tedio=fopen("CIniciales.dat","r")) != NULL) /* file exists */
	fprintf(stderr,"CIniciales existia\n");
	
   else {
     fprintf(stderr,"Guardo transitorio  de %d días\n", TsaveCI);
     tedio=fopen("CIniciales.dat","w");
     for (i=0;i<LIMfilas; i++) {
        for(j=0;j<LIMcolumnas;j++)
	  for(k=0; k< POBLACIONES; k++)
	     fprintf(tedio,"%ld ", p[i][j][k]);
        fprintf(tedio,"\n");
       };
  
     fclose(tedio);
   };
}

void nuevasemilla(long idum) {
   FILE *tedio;
   tedio=fopen("semilla.inp","w");
   fprintf(tedio,"%ld\n", idum);
   fclose(tedio);
}



/* Initial Conditions of the populations */
void pobla_iniciales(int *Transitorio, int *tiempo_inicial) {
  FILE *tedio;
  int i,j,k;
  
  if((tedio=fopen("CIniciales.dat","r")) == NULL) {
    fprintf(stderr,"Transitorio largo\n");
    *tiempo_inicial= TsaveCI;
   /* initialize all in zero less the susceptible mosquitoes and humans susceptibles*/
    for (i=0;i<LIMfilas; i++)
      for(j=0;j < LIMcolumnas;j++) {
        pobla[i][j][0]=(long)rint(C.AmbientCapacity[i][j]);
        pobla[i][j][1]=0;
        pobla[i][j][2]=(long)rint(C.gente[i][j]);
        pobla[i][j][3]=0;
        pobla[i][j][4]=0;
        
      }; 
    /* C.Inf0 infected person in the location C.centro (the value of C.Inf0 is given in CONSTANTES)*/
    pobla[C.centro][C.centro][2]-=C.Inf0;
    pobla[C.centro][C.centro][3]=C.Inf0;
    }
  else { /*reading the initial population from the file*/
    *tiempo_inicial=0.;
    *Transitorio= *Transitorio+TsaveCI;
    for (i=0;i<LIMfilas; i++){
      for(j=0;j<LIMcolumnas;j++){
	for(k=0; k< POBLACIONES; k++)
          fscanf(tedio,"%ld", pobla[i][j]+k);
      }
    };   
  
  }; /*reading the data of tedio */

        
  if((tedio != NULL))
    fclose(tedio); /* closing tedio */
}/* returning with the initial populations */




