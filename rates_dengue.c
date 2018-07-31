/* Here we can found the rates that are not constant in time */
/* also the V(N) spatial distribution */

float rateVuelo(float mosquitos){
   float vuelo;
   //vuelo=0;
   vuelo=0.059;
   if(mosquitos<1800.)
      vuelo+=0.342822071*exp(-0.003777299*mosquitos);
   return vuelo;
}

float rateInfection_ordenCero(long humanosTotales, long humanosInfectados) {
    float rateOrden0;
    
    if(humanosTotales>0)
       rateOrden0=C.bitr*C.trans_hm*humanosInfectados/humanosTotales;
    else
       rateOrden0=0.0;

    return rateOrden0;
}

float rateInfection_PrimerOrden(long humanosTotales, long terminoExtra) {
    float rateOrden1;
    
    if(humanosTotales>0)
       rateOrden1=C.bitr*C.trans_hm*terminoExtra/(humanosTotales*humanosTotales);
    else
       rateOrden1=0.0;

    return rateOrden1;
}

/* V(N) spatial distribution */
   float NumberOfVectors(float N, float Vi){
     float aux;

     /*SIGMOID*/ /*float b,c;
     b=-200.0; 
     c=10000.0;
     float aux_bis;
     aux_bis=sqrt(b*b+c);
     aux=Vi*(aux_bis*(N+b)/sqrt((N+b)*(N+b)+c)-b)/(aux_bis-b);*/
     
     /* lineal */
     aux=Vi*N;
     return aux;
  }

/* Human Distribution */
float DistribucionDeHumanos(int i, int j, int centro){ 
  float gente;
  /* Distribucion Constante*/
  /*gente=1.;*/

  /* You can define the function that you want, for example ...*/
  /*Exponential distribution: increases with C.lambda<0 and decreases with C.lambda>0*/ 
  /* see CONSTANTES to C.lambda value*/
  float aux;
  aux=sqrt((i-centro)*(i-centro)+(j-centro)*(j-centro));
  gente = exp(-C.lambda*aux);
      
    /* other funcions */

  /* up and down (cosine) function whose aplitud decreases (o increases) with the distance to the center depending on lambda sign*/
  /*aux=sqrt((i-centro)*(i-centro)+(j-centro)*(j-centro));
  gente=exp(-C.lambda*aux)*(cos(3.1416*aux/(2*2.236))+2.);*/
  
 /*Function "cuadrado loco"*/
/*float d;
  aux=sqrt((i-centro)*(i-centro)+(j-centro)*(j-centro));
  if(aux<5.0)
     gente=0.5;
  if((aux > 4.5) & (aux < 7.7))
     gente=0.1;
  if((aux > 7.7) & (aux < 10.1))
     gente=0.35;
  if(aux>9.9)
     gente=0.05;*/
  
  return gente;
}


/* Here I give values to the variables C.AmbientCapacity, C.gente and C.vuelo*/

  void DistribucionesEspaciales(int centro,float Vi){
      int i, j;
      float aux,N0,gente[LIMfilas][LIMcolumnas];

      aux=0.0;

     for(i=0;i<LIMfilas;i++)
        for(j=0;j<LIMcolumnas;j++){
           gente[i][j]= DistribucionDeHumanos(i, j, centro);
           aux += gente[i][j];
        }
     
     N0=C.Nh*1./aux;
     
     aux<-0.;
     for(i=0;i<LIMfilas;i++)
        for(j=0;j<LIMcolumnas;j++){
           C.gente[i][j]= N0*gente[i][j];
           aux=NumberOfVectors(N0*gente[i][j], Vi);
           C.AmbientCapacity[i][j]=aux;
           C.vuelo[i][j]=rateVuelo(aux);
           
        }
    
      
   }

/* here rates that depend on temperature */

/* gonotrophic-cycle coefficient */
float ovi(float grados) {
  float load;
  load=(0.04291*exp(((grados)-4.7511)/10.580));
  return load;
};

/* mortality mosquito rate coefficient (not very ok, pay attention)*/
float mortality_mosquito(float grados){
  float load;
  load=(0.08326+exp(-grados/6.));
  return load;
}

/* birth mosquito rate coefficient (not very ok, pay attention)*/
float birth_mosquito(float grados){
  float load;
  load=0.003*grados;
  return load;
}

void coeficientes(float Temp) /* Produces the coefficients that varying per day */
{
    C.mm= mortality_mosquito(Temp);
    C.nm=birth_mosquito(Temp);
    //C.bitr=ovi(Temp);
}

