/*ANTONIO NORELLI PROF. RICCI TERSENGHI*/
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#define W 2*M_PI


/*STRUCT***************************************************/
typedef struct{
	double x, y, vx, vy, dt, mb;
	}crd;
	
typedef struct{
	double u, k, j;
	}ukj;
	
	
/*FUNZIONI************************************************/
crd rkstep(crd);
inline double ax(crd);
inline double ay(crd);
inline ukj pcj(crd);
crd ci(char *);
void errj(crd);

int main(){
	char nfile[10]="ANpton.txt", pto;
	crd p0=ci(&pto);
	crd tmp=p0;
	nfile[5]=pto;
	ukj e;
	double tist=0, ji=pcj(p0).j, dbef2p0=0, dbef1p0=0, dtmpp0=0, per=0, tmax;
	FILE *fp;
	printf("Digitare il tempo totale di integrazione\n");
	scanf("%lf", &tmax);
	int maxfile=1, i=0;
	if(tmax>10000*p0.dt){
		maxfile=(tmax/p0.dt)/10000;
	}
	if((fp=fopen(nfile,"w")) == NULL){
		printf("errore nell'apertura del file");
		}
	fprintf(fp,"#File di dati del punto %c dell'esonero\n#ma=%.4lf mb=%.4lf dt=%lf tmax=%d \n#t........x........y........vx........vy........u.........k.........u+k.......j.........(j(t)-j0)/j0...\n", pto, (1-p0.mb), p0.mb, p0.dt, tmax);
	while(tist<tmax){
		e=pcj(tmp);
		/*CONDIZIONE PER CONTENERE LE DIMENSIONI DEL FILE*/
		if((int)(tist/p0.dt)%maxfile==0){
			fprintf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf %.10lf %.15lf\n", tist, tmp.x, tmp.y, tmp.vx, tmp.vy, e.u, e.k, e.u+e.k, e.j, (e.j-ji)/ji);
		}
		/*CALCOLO DEL PERIODO*/
			dbef2p0=dbef1p0;
			dbef1p0=dtmpp0;
			dtmpp0=sqrt((p0.x-tmp.x)*(p0.x-tmp.x)+(p0.y-tmp.y)*(p0.y-tmp.y));
			if((per==0)&&(dbef1p0<dbef2p0)&&(dbef1p0<dtmpp0)&&(dbef1p0<0.01)){
				per=tist-p0.dt;
				printf("il tempo di ritorno vicino alle condizioni iniziali o il periodo e' %lf\n",per);
				if((pto==56)&&(i==0)){
					i++;
					per=0;
					printf("quello appena stampato e' il primo tempo di ritorno\nvicino le condizioni iniziali,\nil prossimo e' il periodo\n");
				}
			}
		/*PASSO ROUNGE-KUTTA2*/
		tmp=rkstep(tmp);
		tist+=p0.dt;
	}
	fclose(fp);
	if(pto==52){
		errj(p0);
		printf("E' stato creato anche un file contenente gli errori su j in funzione di dt");
	}
}

/*DOPPIO PASSO RUNGE-KUTTA 2*************************************/
crd rkstep(crd k){
	crd plus;
	plus.x=k.x+0.5*k.vx*k.dt;
	plus.y=k.y+0.5*k.vy*k.dt;
	plus.vx=k.vx+0.5*ax(k)*k.dt;
	plus.vy=k.vy+0.5*ay(k)*k.dt;
	plus.dt=k.dt;
	plus.mb=k.mb;

	k.x=k.x+plus.vx*k.dt;
	k.y=k.y+plus.vy*k.dt;
	k.vx=k.vx+ax(plus)*k.dt;
	k.vy=k.vy+ay(plus)*k.dt;
	return k;
}

/*ACCELERAZIONE X E Y********************************************/
inline double ax(crd k){
	return -4*M_PI*M_PI*(((1-k.mb)*(k.x+k.mb))/(pow((-k.mb-k.x)*(-k.mb-k.x)+k.y*k.y,1.5))+k.mb*(k.x-(1-k.mb))/(pow(((1-k.mb)-k.x)*((1-k.mb)-k.x)+k.y*k.y,1.5)))+2*W*k.vy+W*W*k.x;
}
inline double ay(crd k){
	return -4*M_PI*M_PI*k.y*((1-k.mb)/(pow((-k.mb-k.x)*(-k.mb-k.x)+k.y*k.y,1.5))+k.mb/(pow(((1-k.mb)-k.x)*((1-k.mb)-k.x)+k.y*k.y,1.5)))-2*W*k.vx+W*W*k.y;
}

/*ENERGIE U, K, QUANTITA J***************************************/
inline ukj pcj(crd k){
	ukj a;
	a.u=-4*M_PI*M_PI*((1-k.mb)/sqrt((-k.mb-k.x)*(-k.mb-k.x)+k.y*k.y)+k.mb/sqrt(((1-k.mb)-k.x)*((1-k.mb)-k.x)+k.y*k.y));
	a.k=0.5*(k.vx*k.vx+k.vy*k.vy);
	a.j=2*a.u+2*a.k-W*W*(k.x*k.x+k.y*k.y);
	return a;
}

/*ERRORE SU J*****************************************************/
void errj(crd p0){
	p0.dt=0.000001;
	crd k=p0;
	double tistej=0, ji=pcj(k).j, ej;
	FILE *fp1;
	if((fp1=fopen("errj(dt).txt","w")) == NULL){
		printf("errore nell'apertura del file");
		}
	fprintf(fp1,"#(j(t)-j0)/j0 al variare di dt per t=1.5\n#dt.......ej.......logdt.....logej\n");
	while(p0.dt<0.001){
		k=p0;
		tistej=0;
		while(tistej<=1.5){
			k=rkstep(k);
			tistej+=k.dt;
		}
		ej=fabs((pcj(k).j-ji)/ji);
		fprintf(fp1,"%lf %.20lf %lf %lf\n", k.dt, ej, log10(k.dt), log10(ej));
		p0.dt*=5;
	}
	fclose(fp1);
}

/*SCELTA CI*******************************************************/
crd ci(char *a){
	char n=0;
	crd k;
	while((n<48)||(n>56)){
		printf("Digitare il punto dell'esonero desiderato per settare le condizioni iniziali,\nper i punti facoltativi utilizzare 5,6 o 7,\nper il secondo punto dell'ultimo facoltativo digitare 8\nper settare c.i. a piacere digitare 0\n");
		scanf("%c", &n);
	}
	*a=n;
	printf("Digitare il dt (si consiglia 0.000005)\n");
	scanf("%lf", &k.dt);
	if(n==49){
		k.x=0.895;
		k.y=0;
		k.vx=0;
		k.vy=1.018;
		k.mb=0.;
	}else if(n==50){
		printf("Digitare mb (0.0001 esonero, 0.001 Giove, 0.000043 Urano, 0.000003 Terra\n");
		scanf("%lf", &k.mb);
		k.x=0.895;
		k.y=0;
		k.vx=0;
		k.vy=1.018;
	}else if((n==51)||(n==52)){
		k.x=-0.995;
		k.y=0;
		k.vx=0;
		k.vy=0.6;
		k.mb=0.0001;
	}else if(n==53){
		k.x=0.499;
		k.y=0.79;
		k.vx=-0.71421;
		k.vy=0.438729;
		k.mb=0.0001;
	}else if(n==54){
		printf("Digitare mb\n");
		scanf("%lf", &k.mb);
		k.x=0.5-k.mb;
		k.y=sqrt(3./4.);
		k.vx=0;
		k.vy=0.01;
	}else if(n==55){
		k.x=-1.03552;
		k.y=0;
		k.vx=0;
		k.vy=0.331842;
		k.mb=0.0001;
	}else if(n==56){
		k.x=1.06691;
		k.y=0.435202;
		k.vx=0.370381;
		k.vy=-1.79883;
		k.mb=0.0001;
	}else{
		printf("Digitare x(0)\n");
		scanf("%lf", &k.x);
		printf("Digitare y(0)\n");
		scanf("%lf", &k.y);
		printf("Digitare vx(0)\n");
		scanf("%lf", &k.vx);
		printf("Digitare vy(0)\n");
		scanf("%lf", &k.vy);
		printf("Digitare mb\n");
		scanf("%lf", &k.mb);
	}
	return k;
}
	



