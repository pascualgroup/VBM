
/* the Calculation the mean values ​​of the populations with border correction is made (include auxiliary functions). */

void	Gs(double Plus,double Minus,double *G0,double *G1,double *G2)
{/* This function calculates the probabilities to exced by 0, 1 and 2 the population border*/

double z, y, yy;

if(Plus == 0){ *G0 = exp(-Minus);
	      *G1 = *G0*Minus;
	      *G2 = *G1*Minus/2.;
	      return;
	    }/* eventos que suman con tasa cero */

if(Minus == 0){*G0 = exp(-Plus);
		*G1=*G2=0.;
		return;
		} /* eventos que restan con tasa cero */

/* general case */
z=2.  * sqrt(Plus*Minus);
y= z/3.75;
yy= y*y;

if (y < 1) {
*G0=1.0+yy*(3.5156229+yy*(3.0899424+yy*(1.2067492
+yy*(0.2659732+yy*(0.360768e-1+yy*0.45813e-2)))));
*G0= *G0 * exp(-(Plus+Minus));
*G1= 2.*Minus*(0.5+yy*(0.87890594+yy*(0.51498869+yy*(0.15084934
+yy*(0.2658733e-1+yy*(0.301532e-2+yy*0.32411e-3))))));
*G1= *G1 * exp(-(Plus+Minus));
}
else {
y=1/y;
*G0=(exp(z-Plus-Minus)/sqrt(z))*(0.39894228+y*(0.1328592e-1
+y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2
+y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
+y*0.392377e-2))))))));
*G1=0.2282967e-1+y*(-0.2895312e-1+y*(0.1787654e-1
-y*0.420059e-2));
*G1=0.39894228+y*(-0.3988024e-1+y*(-0.362018e-2
+y*(0.163801e-2+y*(-0.1031555e-1+y* *G1))));
*G1 =  sqrt(Minus/Plus)* *G1 * (exp(z-(Plus+Minus))/Plus);
}

*G2 = (- *G1 + Minus * *G0)/Plus;
return;
}

void fix(long n, double Plus, double Minus, double *Norm, double *correc,
double *G0, double *G1, double *G2)
/*This function calculates of the correction by border to the mean value and returns Norma and correc* */ 

{
if((n < 0) ||  (n > 2)) {*Norm=1., *correc=0.; return;}
if ((Plus==0) && (Minus==0))
	{*correc=0.; *Norm=1; *G0=1.; *G1=*G2=0.; /* Evitar 0/0 */
	}
else
	{
Gs(Plus,Minus,G0,G1,G2);
switch (n){
     case 0: {*Norm = 1.-*G0-*G1-*G2;
	if(*Norm <= 0.)fprintf(stderr, "Case 0 Norma negativa  %g %g %g %g %g %g\n", *Norm,*G0, *G1,*G2, Plus, Minus);
	*correc = *G1+2* *G2;
        return; break;};
     case 1: {*Norm = 1- *G1- *G2;
	if(*Norm <= 0.)fprintf(stderr, "Case 1 Norma negativa %g\n", *Norm);
	*correc = *G2;
	return; break;};
     case 2: {*Norm = 1- *G2;
	if(*Norm <= 0.)fprintf(stderr, "Case 2 Norma negativa %g\n", *Norm);
	*correc = 0;
	   }
	  }
	}
}

/* This routine computes the mean values of the populations that are necessary for the rates derivate.*/

void VM(double vm[][LIMcolumnas][POBLACIONES],double Tnl[][LIMcolumnas][EVENTOS])
/* Lambdas, pobla y nevents variables are global. Lbd are the Lambdas */
{
int i,j,k,im,ip,jm,jp,e;
long n;
double LbdP, LbdM, G0, G1, G2, Norm, correc, correc2, hold;

for (i=0 ; i< LIMfilas; i++)
	for (j=0; j < LIMcolumnas; j++)
		{
                 for(e=0; e<EVENTOS;e++)
                      Tnl[i][j][e]=0.0;
/* fix returns the correction by border of the population */

/* fligth indices */
        im=i-1; jm=j-1; ip=i+1;jp=j+1;
        if( (im < 0) ) im=0;
	if( (jm < 0) ) jm=0;
	if( (jp >= LIMcolumnas) ) jp=LIMcolumnas-1;
	if( (ip >= LIMfilas) ) ip=LIMfilas-1;

/* susceptibles mosquitoes*/
	n=pobla[i][j][0];
	LbdM=Lbd[i][j][0]+Lbd[i][j][1]+Lbd[i][j][10];
        LbdP=Lbd[i][j][3]+
          (Lbd[ip][j][10] + Lbd[i][jm][10]+ Lbd[im][j][10] + Lbd[i][jp][10]) *C.pvuel[1]+
	  (Lbd[ip][jp][10]+Lbd[im][jm][10]+ Lbd[ip][jm][10]+ Lbd[im][jp][10])*C.pvuel[0];

        fix(n, LbdP, LbdM, &Norm, &correc, &G0, &G1, &G2);
	vm[i][j][0]= n+(LbdP-LbdM + correc)/Norm;
        Tnl[i][j][0]=LbdP+LbdM-Lbd[i][j][1];

/* Infected mosquitoes */
	n=pobla[i][j][1];
        LbdM=Lbd[i][j][2]+Lbd[i][j][11];
	LbdP=Lbd[i][j][1]+
          (Lbd[ip][j][11] + Lbd[i][jm][11]+ Lbd[im][j][11] + Lbd[i][jp][11]) *C.pvuel[1]+
	  (Lbd[ip][jp][11]+Lbd[im][jm][11]+ Lbd[ip][jm][11]+ Lbd[im][jp][11])*C.pvuel[0]; 
	
	fix(n, LbdP, LbdM, &Norm, &correc, &G0, &G1, &G2);
	vm[i][j][1]= n+(LbdP-LbdM + correc)/Norm;
        Tnl[i][j][2]=LbdP+LbdM-Lbd[i][j][1];

/* Suceptibles humans */
	n=pobla[i][j][2];
	LbdP=Lbd[i][j][9];
	LbdM=Lbd[i][j][4]+Lbd[i][j][5];
	fix(n, LbdP, LbdM, &Norm, &correc, &G0, &G1, &G2);
	vm[i][j][2]= n+(LbdP-LbdM + correc)/Norm;
        
        Tnl[i][j][1]-=pobla[i][j][3]*(LbdP-LbdM + correc)/Norm;
        Tnl[i][j][5]+=(pobla[i][j][3]+pobla[i][j][4])*(LbdP-LbdM + correc)/Norm;
       
        

/* Infected humans */
	n=pobla[i][j][3];
	LbdP=Lbd[i][j][5];
	LbdM=Lbd[i][j][6]+Lbd[i][j][7];
	fix(n, LbdP, LbdM, &Norm, &correc, &G0, &G1, &G2);
	vm[i][j][3]= n+(LbdP-LbdM + correc)/Norm;
     
        Tnl[i][j][1]+=(pobla[i][j][2]+pobla[i][j][4])*(LbdP-LbdM + correc)/Norm;
        Tnl[i][j][5]-=pobla[i][j][2]*(LbdP-LbdM + correc)/Norm;

/* Recovered humans */
	n=pobla[i][j][4];
	LbdP=Lbd[i][j][7];
	LbdM=Lbd[i][j][8];
	fix(n, LbdP, LbdM, &Norm, &correc, &G0, &G1, &G2);
	vm[i][j][4]= n+(LbdP-LbdM + correc)/Norm;
        
        Tnl[i][j][1]-=pobla[i][j][3]*(LbdP-LbdM + correc)/Norm;
        Tnl[i][j][5]-=pobla[i][j][2]*(LbdP-LbdM + correc)/Norm;
}; // end of i j
} 
