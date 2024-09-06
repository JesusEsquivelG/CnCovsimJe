/*
Compile 
gcc SIQR_COV_PAPER_RSTA_SEPT24_fig4.c -lm -o SIQR_COV_PAPER_RSTA_SEPT24_fig4
Run 
nohup ./SIQR_COV_PAPER_RSTA_SEPT24_fig4 10000 1000 2 0.3 0 1 0 0.2213 5 5 10 1 100 &  
nohup ./SIQR_COV_PAPER_RSTA_SEPT24_fig4 10000 1000 2 0.7 0 1 0 0.2213 5 5 10 1 100 & 
nohup ./SIQR_COV_PAPER_RSTA_SEPT24_fig4 10000 1000 2 1 0 1 0 0.2213 5 5 10 1 100 &
nohup ./SIQR_COV_PAPER_RSTA_SEPT24_fig4 10000 1000 2 1 2 1 0 0.2213 5 5 10 1 100 &
*/

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <time.h>
#include <float.h>

int getAge ( void );
void spreadDiseaseFinal (char nameDir[], char nameFiles[], unsigned int numNodes, unsigned int *LISTA_VECINOS[], unsigned int *ARRAY_DEGREE, unsigned int numSimulations, unsigned int spreadSteps, unsigned int I0, float SI,  int IQ, int QR, int IR, float S[], float I[], float Q[], float R[], int AgeNod[], float *ageNodInfectedRate);
int array_search_1(unsigned int array[], int size, int dato);
void random_uniform(unsigned int M,unsigned long N,unsigned long *num);  
float random_uniform_f() ;
void calculateSimplicesFinal (unsigned int *ARRAY_DEGREE, unsigned int *LISTA_VECINOS[], unsigned int numNodes, char nameDir[], char nameFiles[]);
void meanDegreeCalculation (unsigned int connectedNodes, float meanDegree[], unsigned int *ARRAY_DEGREE, unsigned int simulations);
int getNumLinks (unsigned int idNode, int AgeNod[], unsigned int connectedNodes, unsigned int *ARRAY_DEGREE, float P, unsigned int initialLinks);
void preferentialLinking (unsigned int idNode, int AgeNod[], float P, unsigned int initialLinks,  unsigned int connectedNodes, unsigned int Ai, unsigned int div_Ai, unsigned int *ARRAY_DEGREE, unsigned int *LISTA_VECINOS[]);
void getPk (char nameDir[], char nameFiles[], unsigned int numNodes, unsigned int *ARRAY_DEGREE, unsigned int *LISTA_VECINOS[]);
void connectNodes (unsigned int *LISTA_VECINOS[], unsigned int *ARRAY_DEGREE, unsigned int sNode, unsigned int dNode);

double totalIntfected=0.0;
float beta[2000];//To calculate the CDF of beta = psi/<k> of the paper up to 2000 simulations
unsigned int debugComments=0, simNumber=0;

int main(int argc, char *argv[]) {
	
unsigned int Ai, div_Ai, initialLinks, numNodes, numSimulations, z=0, numSteps, spreadSteps, I0, quarantineRate, recoveryRateAfterQuarantine, recoveryRateWithoutQuar;
float P, infectionRate;
char numSimulationsString[10]="", numNodesString[10]="";
char nameDir[620]="ResCovAgeMexJun24_";
char nameFiles[620]="_";
float meanDegree[100000];
for (z=0;z<1000;z++) meanDegree[z]=0.0;
int AgeNod[10000];//Age group of the node
float ageNodRate[22];
float ageNodInfectedRate[22];

if (argc != 14)  {
	printf("-------> ERROR! INCOMPLETE DATA %d != 14\n", argc);
} else {
	numNodes = (int) strtol(argv[1], NULL, 10);         
	numSimulations = (int) strtol(argv[2], NULL, 10);
	initialLinks = (int) strtol(argv[3], NULL, 10);
	P= 1.0*(float) strtof(argv[4], NULL);
	Ai = (int) strtol(argv[5], NULL, 10);
	div_Ai = (int) strtol(argv[6], NULL, 10);
	debugComments = (int) strtol(argv[7], NULL, 10);
	infectionRate = 1.0*((float) strtof(argv[8], NULL));
	quarantineRate = ((int) strtol(argv[9], NULL, 10));
	recoveryRateAfterQuarantine = ((int) strtol(argv[10], NULL, 10));
	recoveryRateWithoutQuar = ((int) strtol(argv[11], NULL, 10));
	I0 = ((int)	  strtol(argv[12], NULL, 10));
	spreadSteps = ((int) strtol(argv[13], NULL, 10));    
	if (debugComments == 1) {
		printf("SALIDA A TERMINAL PARA DEPURACION: ACTIVADA\n");	
	} else {
		printf("SALIDA A TERMINAL PARA DEPURACION: DESACTIVADA\n");	
	}

	strcat(numNodesString, argv[1]);
	strcat(numSimulationsString, argv[2]);
	//STRINGS FOR FILE NAMES
	unsigned int a=0;
	char dir[500] = "mkdir ";
	strcat(nameDir, argv[1]);strcat(nameDir, "n_");
	strcat(nameDir, "copyModel_");	
	strcat(nameDir, argv[1]);strcat(nameDir, "n_m_");
	strcat(nameDir, argv[3]);strcat(nameDir, "_P_");
	strcat(nameFiles, "copyModel_");
	strcat(nameFiles, argv[1]);strcat(nameFiles, "n_m_");
	strcat(nameFiles, argv[3]);strcat(nameFiles, "_P_");
	strcat(nameDir, argv[4]);strcat(nameDir, "_Ai_");strcat(nameDir, argv[5]);strcat(nameDir, "_");strcat(nameDir, argv[6]);
	strcat(nameFiles, argv[4]);strcat(nameFiles, "_Ai_");strcat(nameFiles, argv[5]);strcat(nameFiles, "_");strcat(nameFiles, argv[6]);
	strcat(nameDir, "_psi_");strcat(nameDir, argv[8]);
	strcat(nameFiles, "_psi_");strcat(nameFiles, argv[8]);
	strcat(nameDir, "_sigma_");strcat(nameDir, argv[9]);
	strcat(nameFiles, "sigma");strcat(nameFiles, argv[9]);
	strcat(nameDir, "_epsilon_");strcat(nameDir, argv[10]);
	strcat(nameFiles, "_epsilon_");strcat(nameFiles, argv[10]);
	strcat(nameDir, "_eta_");strcat(nameDir, argv[11]);
	strcat(nameFiles, "_eta_");strcat(nameFiles, argv[11]);
	strcat(nameDir, "_I0_");strcat(nameDir, argv[12]);
	strcat(nameFiles, "_I0_");strcat(nameFiles, argv[12]);
	strcat(nameDir, "_spreadSteps_");strcat(nameDir, argv[13]);
	strcat(nameFiles, "_spreadSteps_");strcat(nameFiles, argv[13]);		
	strcat(dir, nameDir);
	system(dir);//se crea el directorio de esta simulacion


	unsigned int *ARRAY_DEGREE;
	unsigned int *ARRAY_STRENGTH;
	unsigned int *LISTA_VECINOS[numNodes];
	unsigned int *LISTA_VECINOS_COMMUNITY[numNodes];

	unsigned int n=0;
	float S[spreadSteps];
	float I[spreadSteps];
	float Q[spreadSteps];
	float R[spreadSteps];

	for (z=0;z<spreadSteps;z++) {
		S[z]=0.0;
		I[z]=0.0;
		Q[z]=0.0;
		R[z]=0.0;
	}

/*****************************************************/
struct timespec ts1; 
clock_gettime( CLOCK_REALTIME, &ts1 );
srandom(ts1.tv_nsec);
/******************************************************************************************/
	for(z=0;z<numSimulations;z++) {
		printf("\n*******************SIMULATION  %d***************************\n",z);
		//CICLO PARA RESERVAR LA MEMORIA DINAMICA PARA EL ARREGLO DE LOS VECINOS
		for (n=0;n<numNodes;n++) {LISTA_VECINOS[n] = calloc(2, sizeof(int));if (LISTA_VECINOS[n]  == NULL) {printf("Error al intentar reservar memoria\n");}}		
		ARRAY_DEGREE = calloc(numNodes, sizeof(int));if (ARRAY_DEGREE == NULL) {printf("Error al intentar reservar memoria\n");}
		ARRAY_STRENGTH = calloc(numNodes, sizeof(int));if (ARRAY_STRENGTH == NULL) {printf("Error al intentar reservar memoria\n");}	
		for(n=0;n<numNodes;n++) {
		 	ARRAY_DEGREE[n]=0;
		}    
		//models with addition of nodes and links at each time step
		//Three fully connected nodes at start
		unsigned int connectedNodes = 3;
		AgeNod[0]=getAge();		
		AgeNod[1]=getAge();
		AgeNod[2]=getAge();
		ageNodRate[AgeNod[0]]+=1.0;
		ageNodRate[AgeNod[1]]+=1.0;
		ageNodRate[AgeNod[2]]+=1.0;
		connectNodes (LISTA_VECINOS, ARRAY_DEGREE, 1, 0);
		connectNodes (LISTA_VECINOS, ARRAY_DEGREE, 2, 0);
		connectNodes (LISTA_VECINOS, ARRAY_DEGREE, 2, 1);					
		//LOOP FOR ADD NODES
		int k=0;
		for(k=connectedNodes;k<numNodes;k++) {
			meanDegreeCalculation (connectedNodes, meanDegree, ARRAY_DEGREE, numSimulations);	
			if (debugComments) printf("\n***************************************************>  NODE %d WILL BE ADDED\n",k);
			AgeNod[k]=getAge();
			ageNodRate[AgeNod[k]]+=1.0;
			preferentialLinking (k, AgeNod, P, initialLinks, connectedNodes, Ai, div_Ai, ARRAY_DEGREE, LISTA_VECINOS);
			connectedNodes++;
		}
		getPk(nameDir, nameFiles, numNodes, ARRAY_DEGREE, LISTA_VECINOS);
		double kPromSim = 0;	
		int q=0;
		for (q=0;q<numNodes;q++) kPromSim+=1.0*ARRAY_DEGREE[q]/numNodes;
		float SI_=1.0*infectionRate;
		// float SI_=1.0*infectionRate/kPromSim;
		// beta[z]=SI_;
		spreadDiseaseFinal (nameDir, nameFiles, numNodes, LISTA_VECINOS, ARRAY_DEGREE, numSimulations, spreadSteps, I0, SI_, quarantineRate, recoveryRateAfterQuarantine, recoveryRateWithoutQuar, S, I, Q, R, AgeNod, ageNodInfectedRate);
		printf("\n***************************************************>  beta = %f/%f=%f\n", 1.0*infectionRate, kPromSim,SI_);
		int h=0, v=0;
		//liberamos la memoria
		for (n=0;n<numNodes;n++) {
			free(LISTA_VECINOS[n]);        
		} 	       
	}

	/*************************** betaCDF **************************/
	int ub = 0, vb=0, indiceb=0;
	float betaOr[numSimulations];
	for (vb=0;vb<numSimulations;vb++) {
		printf("Sim %d, beta=%f \n", vb, beta[vb]);
		betaOr[vb]=beta[vb];
	}
	float valCb = 0.0;
	float minCb = 1000.0;
	float Beta_redOrdenado[numSimulations];
	float sumaCb = 0.0;
	for (vb=0;vb<numSimulations;vb++) {
		minCb=1000;	
		for (ub=0;ub<numSimulations;ub++) {
			if (beta[ub]<minCb) {
				minCb = beta[ub];
				indiceb=ub;
			}
		}
		Beta_redOrdenado[vb]=minCb;
		beta[indiceb]=10000;
		indiceb=0;	
	}
	for (vb=0;vb<numSimulations;vb++) sumaCb=sumaCb+Beta_redOrdenado[vb];
	char cadCb[1000]="";
	strcat(cadCb,nameDir);	
	strcat(cadCb,"/");	
	strcat(cadCb,"CDF_beta");
	strcat(cadCb,nameFiles);
	strcat(cadCb,".txt"); 
	FILE *ficheroCb;
	ficheroCb = fopen(cadCb, "a+");
	for (vb=0;vb<numSimulations;vb++) {
		float suma = 0.0;
		for (ub=0;ub<=vb;ub++) suma+=1.0*Beta_redOrdenado[ub]/sumaCb;
		fprintf(ficheroCb, "%f\t%f\t%f\n", Beta_redOrdenado[vb], 1.0*suma, betaOr[vb]);
		//printf("%f\t%f\n", C_redOrdenado[v],1.0*suma);
	}
	fclose(ficheroCb); 

	/*************************** P(k) **************************/
	char cad03[400]="";
	strcat(cad03,nameDir);	
	strcat(cad03,"/");	
	strcat(cad03,"dist_grado");
	strcat(cad03,nameFiles);
	strcat(cad03,".txt");
	char enlace_out[870]="tclsh Promedio.tcl ";
	char dat_out[270]="Datos_";
	strcat(enlace_out,cad03);
	strcat(enlace_out," ");
	strcat(enlace_out,numNodesString);
	strcat(enlace_out," ");
	strcat(enlace_out,numSimulationsString);
	strcat(enlace_out," >> ");
	strcat(enlace_out,nameDir);	
	strcat(enlace_out,"/");	
	strcat(dat_out,"dist_grado");
	strcat(dat_out,nameFiles);
	strcat(dat_out,".txt");
	strcat(enlace_out,dat_out);
	printf("**************** %s\n", enlace_out);  
	system(enlace_out);

	//PARA OBTENER LOS PROMEDIOS S,I,Q,R POR ITERACION
	char cadt[850]="";
	strcat(cadt,nameDir);	
	strcat(cadt,"/");	
	strcat(cadt,"S_I_Q_R_");
	strcat(cadt,nameFiles);
	strcat(cadt,".txt");
	printf("S \t I\t Q\t R\t \n");  
	FILE *ficherot;
	ficherot = fopen(cadt, "a+");
	for (a=0;a<spreadSteps;a++) {
	//	printf("%f \t %f\t %f\t %f\t \n", S[a], I[a], Q[a], R[a]);  
				fprintf(ficherot, "%d\t %f \t %f\t %f\t %f\t %f\n", a, 1.0*S[a]/numNodes, 1.0*I[a]/numNodes, 1.0*Q[a]/numNodes, 1.0*R[a]/numNodes, 1.0*(S[a]+I[a]+Q[a]+R[a])/numNodes);
//				fprintf(ficherot, "%d\t %f \t %f\t %f\t %f\t %f\n", a, S[a], I[a], Q[a], R[a], S[a]+I[a]+Q[a]+R[a]);
	}
	fclose(ficherot); 	

	char xticsLabels[23][6] = {"0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54", "55-59", "60-64", "65-69", "70-74", "75-79", "80-84", "85-89", "90-94", "95-99", "100-"};
	//age
	char cadage[850]="";
	strcat(cadage,nameDir);	
	strcat(cadage,"/");	
	strcat(cadage,"AGE_");
	strcat(cadage,nameFiles);
	strcat(cadage,".txt");
	FILE *ficheroage;
	ficheroage = fopen(cadage, "a+");
	for (a=0;a<22;a++) {
		fprintf(ficheroage, "%d\t %f\t %s\n", a, 1.0*ageNodRate[a]/(numSimulations*numNodes), xticsLabels[a]);
	}
	fclose(ficheroage); 


	//infected by age
	char cadageI[850]="";
	strcat(cadageI,nameDir);	
	strcat(cadageI,"/");	
	strcat(cadageI,"AGE_INFECTED_");
	strcat(cadageI,nameFiles);
	strcat(cadageI,".txt");
	FILE *ficheroageI;
	ficheroageI = fopen(cadageI, "a+");
	for (a=0;a<22;a++) {
		fprintf(ficheroageI, "%d\t %f\t %s\n", a, 1.0*ageNodInfectedRate[a]/(numSimulations*numNodes), xticsLabels[a]);
	}
	fclose(ficheroageI); 	

	//<k>
	char cadgp[850]="";
	strcat(cadgp,nameDir);	
	strcat(cadgp,"/");	
	strcat(cadgp,"Kprom_");
	strcat(cadgp,nameFiles);
	strcat(cadgp,".txt");
	FILE *ficherogp;
	ficherogp = fopen(cadgp, "a+");
	for (a=0;a<100000;a++) {
		fprintf(ficherogp, "%d\t", a*1);
		fprintf(ficherogp, "%f \t", meanDegree[a]);
		fprintf(ficherogp, "\n");		
	}
	fclose(ficherogp);

	//Plot File
	char cadplot[850]="";
	strcat(cadplot,nameDir);	
	strcat(cadplot,"/");	
	strcat(cadplot,"Plot_");
	strcat(cadplot,nameFiles);
	strcat(cadplot,".txt");
	FILE *ficheroPlot;
	ficheroPlot = fopen(cadplot, "a+");
		fprintf(ficheroPlot, "set grid\n");
		fprintf(ficheroPlot, "set title \" ");
		fprintf(ficheroPlot, "Model Copy P=%f, m=%d, Ai=%f\"\n", P, initialLinks, 1.0*(Ai/div_Ai));

		fprintf(ficheroPlot, "set multiplot layout 2,2\n");
		fprintf(ficheroPlot, "set logscale\n");
		fprintf(ficheroPlot, "plot ");
		fprintf(ficheroPlot, " \"%s/%s\" using 1:2 title \"Degree distribution P(k)\" with lp lt 1 lw 1 pt 6 lc \"black\" \n", nameDir, dat_out);

		// fprintf(ficheroPlot, "set title \" Community Size distribution P(s) \n");
		// fprintf(ficheroPlot, "plot ");
		// fprintf(ficheroPlot, " \"%s/%s\" using 1:2 title \"P(s) \" with lp lt 1 lw 1 pt 7 lc \"black\" \n", nameDir, dat_outSc);

		fprintf(ficheroPlot, "unset logscale\n");	
		fprintf(ficheroPlot, "set xtics rotate by -65\n");
		fprintf(ficheroPlot, "set title \"AGE INFECTED (TOTAL INFECTED=%f)\n", totalIntfected);
		fprintf(ficheroPlot, "plot ");
		fprintf(ficheroPlot, " \"%s\" using 1:2:xticlabels(3) title \"AGE INFECTED\" with lp lt 1 lw 2 pt 4 lc \"black\" \n", cadageI);

		// fprintf(ficheroPlot, "set title \" <k> \n");
		// fprintf(ficheroPlot, "plot [100:%d]", numNodes-100);
		// fprintf(ficheroPlot, " \"%s\" using 1:2 title \"<k>\" with lp lt 1 lw 2 pt 4 lc \"black\" \n", cadgp);
		fprintf(ficheroPlot, "set title \"AGE \n");
		fprintf(ficheroPlot, "plot ");
		fprintf(ficheroPlot, " \"%s\" using 1:2:xticlabels(3) title \"AGE\" with lp lt 1 lw 2 pt 4 lc \"black\" \n", cadage);


		fprintf(ficheroPlot, "set title \" Health dynamic ");
		fprintf(ficheroPlot, "psi: %f, ", infectionRate);
		fprintf(ficheroPlot, "sigma: %d days, ", quarantineRate);
		fprintf(ficheroPlot, "epsilon %d days, ", recoveryRateAfterQuarantine);
		fprintf(ficheroPlot, "eta: %d days, ", recoveryRateWithoutQuar);
		if (I0>0) fprintf(ficheroPlot, "I0: %d, ", I0);
		fprintf(ficheroPlot, " \"\n");

		fprintf(ficheroPlot, "plot [1:150]");
		fprintf(ficheroPlot, " \"%s\" using 1:2 title \"S\" with lines lw 2 lc \"green\", ", cadt);
		fprintf(ficheroPlot, " \"%s\" using 1:3 title \"I\" with lines lw 2 lc \"red\", ", cadt);
		fprintf(ficheroPlot, " \"%s\" using 1:4 title \"Q\" with lines lw 2 lc rgb \"#FF00FF\", ", cadt);
		fprintf(ficheroPlot, " \"%s\" using 1:5 title \"R\" with lines lw 2 lc \"blue\" \n", cadt);	
		fprintf(ficheroPlot, "\nunset multiplot\n");
	fclose(ficheroPlot);


	char plot[1000]="gnuplot -persist ";
	strcat(plot,cadplot);
	system(plot);

	}
}
/****************************SPREAD DISEASE*******************************/
void spreadDiseaseFinal (char nameDir[], char nameFiles[], unsigned int numNodes, unsigned int *LISTA_VECINOS[], unsigned int *ARRAY_DEGREE, unsigned int numSimulations, unsigned int spreadSteps, unsigned int I0, float infectionR, int quarantineR, int recoveryAfterQ, int recoveryWithoutQ, float S[], float I[], float Q[], float R[], int AgeNod[], float *ageNodInfectedRate) {
	//obtenemos el promedio de grado
    unsigned long random;
    float num_r;
    int k=0,s=0;
    //PARA OBTENER EL STATUS FINAL DE CADA NODO en cada simulacion 
    char cadStatusFinal[850]="";
    strcat(cadStatusFinal,nameDir);	
    strcat(cadStatusFinal,"/");	
    strcat(cadStatusFinal,"HEALTH_FINAL");
    strcat(cadStatusFinal,nameFiles);
    strcat(cadStatusFinal,".txt");
    printf("\n*******************INICIO PROPAGACION DE INFECCION SI=%f***************************\n", infectionR);
    //INICIAMS EL PROCESO DE INFECCION
    //CREAMOS UN ARREGLO PARA CADA ESTADO DE LOS NODOS
    char ESTADO[10000];//[CANTIDAD_DE_NODOS];
    char ESTADO_TEMPORAL[10000];//[CANTIDAD_DE_NODOS];
    unsigned int DAYS_TO_RECOVERY_INFECTED[10000];
	unsigned int DAYS_TO_GETOUT_QUARANTINE[10000];
    int CANTIDAD_INFECTADOS=0,CANTIDAD_SUCEPTIBLES=numNodes,CANTIDAD_EN_CUARENTENA=0,CANTIDAD_REMOVIDOS=0;

    if (debugComments) printf("COMENZARA LA EVOLUCION DE LA INFECCION\n");
    if (debugComments) printf("PROB S->I = %f\n", infectionR);
    if (debugComments) printf("PROB I->Q = %f\n", 1.0/quarantineR);
    if (debugComments) printf("PROB Q->R = %f\n", 1.0/recoveryAfterQ);
    if (debugComments) printf("PROB I->R = %f\n", 1.0/recoveryWithoutQ);

    /*MODELO:
    * A CADA PASO DE TIEMPO:
    * CON PROBABILIDAD SI UN NODO SUCEPTIBLE QUE ESTE CONECTADO A UN N9ODO INFECTADO SE INFECTA.
    * CON PROBABILIDAD SR UN NODO SUCEPTIBLE ES VACUNADO
    * CON PROBABILIDAD IQ UN NODO INFECTADO ES PUESTO EN CUARENTENA
    * CON PROBABILIDAD IR UN NODO INFECTADO SE RECUPERA
    * CON PROBABILIDAD QR UN NODO EN CUARENTENA SE RECUPERA
    * CON PROBABILIDAD RS UN NODO RECUPERADO SE VUELVE SUCEPTIBLE
    */
    if (debugComments) {
        int h=0, v=0;
        //IMPRIMIMOS LA LISTA DE VECINOS DE CADA NODO
        for (h=0;h<numNodes;h++) {
            printf("++++++++++++LOS VECINOS DEL NODO %d SON %d:\n", h, ARRAY_DEGREE[h]);
            for (v=0;v<ARRAY_DEGREE[h];v++) {
                printf("%d,",LISTA_VECINOS[h][v]);				
            }		
            printf("\n");
        }
    }
    int x=0;
    for (x=0;x<numNodes;x++) {
        ESTADO[x]='S';
        ESTADO_TEMPORAL[x]='S';        
        DAYS_TO_RECOVERY_INFECTED[x]=0;
		DAYS_TO_GETOUT_QUARANTINE[x]=0;
    }
	int totalInfected=0;
	num_r=random_uniform_f();
    for (x=0;x<spreadSteps;x++) {
        for(k=0;k<numNodes;k++) {
            if (ESTADO[k]=='I' && DAYS_TO_RECOVERY_INFECTED[k]==0) {
                    ESTADO_TEMPORAL[k]='R';
                    CANTIDAD_REMOVIDOS++;
                    CANTIDAD_INFECTADOS--;
                    if (debugComments) printf("EL NODO INF %d SE RECUPERO\n", k); 
            } else if (ESTADO[k]=='Q' && DAYS_TO_RECOVERY_INFECTED[k]==0) {
                    ESTADO_TEMPORAL[k]='R';
                    CANTIDAD_EN_CUARENTENA--;
                    CANTIDAD_REMOVIDOS++;
                    if (debugComments) printf("EL NODO QUA %d SE RECUPERO\n", k); 
            } else if (ESTADO[k]=='I' && num_r<1.0/quarantineR) {
                    ESTADO_TEMPORAL[k]='Q';
                    CANTIDAD_EN_CUARENTENA++;
                    CANTIDAD_INFECTADOS--;
					DAYS_TO_RECOVERY_INFECTED[k]--;
                    if (debugComments) printf("EL NODO %d SE PUSO EN QUARENTENA\n", k); 
            } else  if (ESTADO[k]=='I') {
				DAYS_TO_RECOVERY_INFECTED[k]--;
				ESTADO_TEMPORAL[k]='I';
				if (debugComments) printf("EL NODO %d TIENE %d DIAS INFECTADO\n", k, 15-DAYS_TO_RECOVERY_INFECTED[k]); 
			}  else  if (ESTADO[k]=='Q') {
				DAYS_TO_RECOVERY_INFECTED[k]--;
				ESTADO_TEMPORAL[k]='Q';
				if (debugComments) printf("EL NODO %d TIENE %d DIAS en cuarenena\n", k, 15-DAYS_TO_RECOVERY_INFECTED[k]); 
			}
        }          
        if (debugComments) printf("ITERACION EPIDEMICA %d\n", x); 
        int ii=0, h=0, v=0, s=0;		
        if (x==0) {
            //SE INFECTAN I0 NODOS AL HAZAR
            while (ii<I0) {			
                random_uniform(0,numNodes-1,&random);	
                if (ESTADO[random]!='I') {// && ARRAY_DEGREE[random]<14) {
                    ESTADO_TEMPORAL[random]='I';
                    ii++;
                    CANTIDAD_SUCEPTIBLES--;
                    CANTIDAD_INFECTADOS++;
					totalInfected++;
					ageNodInfectedRate[AgeNod[random]]+=1.0;					
                    DAYS_TO_RECOVERY_INFECTED[random]=15;
                    if (debugComments) printf("EL NODO %ld HA SIDO INFECTADO INICIALMENTE CON %d DIAS COMO INFECTADO\n", random, DAYS_TO_RECOVERY_INFECTED[random]);
                } else {
					printf("degree = %d\n", ARRAY_DEGREE[random]);
				}
            }
        } else {
			k=0, s=0;
			for(k=0;k<numNodes;k++) {
				//printf("NODO %d = %c", k, ESTADO[k]);
				if (ESTADO[k]=='S') {
					num_r=random_uniform_f();
					int senal=0;
						//RECORREMOS A LOS VECINOS DEL NODO K PARA VER SI ALGUNO ESTA INFECTADO
						for (v=0;v<ARRAY_DEGREE[k];v++) {
							//printf("EL NODO %d ES VECINO DEL NODO %d Y SU HEALTH ES: %c\n",  LISTA_VECINOS[k][v], k, ESTADO[ LISTA_VECINOS[k][v]]);								
							num_r=random_uniform_f();
							if (ESTADO[LISTA_VECINOS[k][v]]=='I' && num_r<infectionR) {																											
							//if (ESTADO[LISTA_VECINOS[k][v]]=='I' && num_r<1.0*covAge[AgeNod[k]]) {																												
								ESTADO_TEMPORAL[k]='I';
								senal=1;
								if (debugComments) printf("EL NODO %d ES INFECTADO POR EL NODO %d\n", k, LISTA_VECINOS[k][v]);
								CANTIDAD_SUCEPTIBLES--;
								CANTIDAD_INFECTADOS++;		
								totalInfected++;
								ageNodInfectedRate[AgeNod[random]]+=1.0;
								DAYS_TO_RECOVERY_INFECTED[k]=recoveryWithoutQ;
								break;						
							}
						}							
					if (senal==0) {
						ESTADO_TEMPORAL[k]='S';
					}
				}
			}              
        }
        k=0;
        for (k=0;k<10000;k++) ESTADO[k]=ESTADO_TEMPORAL[k];
        S[x]+=1.0*(1.0*CANTIDAD_SUCEPTIBLES/numSimulations);///numNodes;
        I[x]+=1.0*(1.0*CANTIDAD_INFECTADOS/numSimulations);///numNodes;
        Q[x]+=1.0*(1.0*CANTIDAD_EN_CUARENTENA/numSimulations);///numNodes;
        R[x]+=1.0*(1.0*CANTIDAD_REMOVIDOS/numSimulations);///numNodes;				
        printf("S=%d\t I=%d\t Q=%d\t R=%d \n", CANTIDAD_SUCEPTIBLES, CANTIDAD_INFECTADOS, CANTIDAD_EN_CUARENTENA, CANTIDAD_REMOVIDOS);
        if (numNodes!=(CANTIDAD_SUCEPTIBLES+CANTIDAD_INFECTADOS+CANTIDAD_EN_CUARENTENA+CANTIDAD_REMOVIDOS)) {
            printf("ERROR!!!!!");
            printf("ERROR!!!!!");
            printf("ERROR!!!!!");
            printf("ERROR!!!!!");
            while(1);
        }
    }	
	printf("TOTALINFECTED=%d\n",totalInfected);
	totalIntfected+=1.0*(1.0*totalInfected/(numNodes*numSimulations));
}
/****************************CONNECT NODES FUNCTION*******************************/
void connectNodes (unsigned int *LISTA_VECINOS[], unsigned int *ARRAY_DEGREE, unsigned int sNode, unsigned int dNode) {
	if (array_search_1(LISTA_VECINOS[sNode], ARRAY_DEGREE[sNode], dNode)==0 && array_search_1(LISTA_VECINOS[dNode], ARRAY_DEGREE[dNode], sNode)==0 && sNode!=dNode) {
		if (debugComments) printf ("EL NODO %d SE CONECTO AL NODO %d\n", sNode, dNode);
		ARRAY_DEGREE[sNode]=ARRAY_DEGREE[sNode]+1;
		void *tmp_ptr3 = realloc(LISTA_VECINOS[sNode], (ARRAY_DEGREE[sNode]+1)*sizeof(int));
			if (tmp_ptr3 == NULL) {
				printf("tomar medidas necesarias\n");
			} else {
			//printf("Reasignación exitosa\n");
				LISTA_VECINOS[sNode] = tmp_ptr3;
			}
		/******************************/  
		LISTA_VECINOS[sNode][ARRAY_DEGREE[sNode]-1]=dNode;


		ARRAY_DEGREE[dNode]=ARRAY_DEGREE[dNode]+1;
		void *tmp_ptr4 = realloc(LISTA_VECINOS[dNode], (ARRAY_DEGREE[dNode]+1)*sizeof(int));
			if (tmp_ptr4 == NULL) {
				printf("tomar medidas necesarias\n");
			} else {
			//printf("Reasignación exitosa\n");
				LISTA_VECINOS[dNode] = tmp_ptr4;
			}
		/******************************/  
		LISTA_VECINOS[dNode][ARRAY_DEGREE[dNode]-1]=sNode;
		//if (debugComments) printf ("EL NODO %d SE CONECTO AL NODO %d\n", dNode, sNode);	
		//ORDENARE LOS VECINOS DEL NODOS 
		unsigned int *LISTA_VEC_AUX1;
		LISTA_VEC_AUX1 = calloc(ARRAY_DEGREE[sNode], sizeof(int));if (LISTA_VEC_AUX1 == NULL) {printf("Error al intentar reservar memoria\n");}		
		int x=0, y=0;
		for (x=0;x<sNode;x++) {
			if (array_search_1(LISTA_VECINOS[sNode], ARRAY_DEGREE[sNode], x)==1) {
				LISTA_VEC_AUX1[y]=x;
				y++;
			}
		}
		for (x=0;x<ARRAY_DEGREE[sNode];x++) LISTA_VECINOS[sNode][x]=LISTA_VEC_AUX1[x];
		
		unsigned int *LISTA_VEC_AUX2;
		LISTA_VEC_AUX2 = calloc(ARRAY_DEGREE[dNode], sizeof(int));if (LISTA_VEC_AUX2 == NULL) {printf("Error al intentar reservar memoria\n");}		
		x=0, y=0;
		for (x=0;x<=sNode;x++) {
			if (array_search_1(LISTA_VECINOS[dNode], ARRAY_DEGREE[dNode], x)==1) {
				LISTA_VEC_AUX2[y]=x;
				y++;
			}
		}
		for (x=0;x<ARRAY_DEGREE[dNode];x++) LISTA_VECINOS[dNode][x]=LISTA_VEC_AUX2[x];
		free(LISTA_VEC_AUX1);
		free(LISTA_VEC_AUX2);
	}
}
/********************************* DEGREE DISTRIBUTION *************************/
void getPk (char nameDir[], char nameFiles[], unsigned int numNodes, unsigned int *ARRAY_DEGREE, unsigned int *LISTA_VECINOS[]) {
	//printf("\n*************************************************Pk");
	char cad[800]="";
	strcat(cad,nameDir);	
	strcat(cad,"/");	
	strcat(cad,"dist_grado");
	strcat(cad,nameFiles);
	strcat(cad,".txt");
	if (debugComments) {
		int h=0, v=0;
		//IMPRIMIMOS LA LISTA DE VECINOS DE CADA NODO
		for (h=0;h<numNodes;h++) {
			printf("---------------LOS VECINOS DEL NODO %d SON:\n", h);
			for (v=0;v<ARRAY_DEGREE[h];v++) {
				printf("%d,",LISTA_VECINOS[h][v]);				
			}		
			printf("\n");
		}
	}	
	//printf("\n************************************************\n");

	if (debugComments)	printf("%s\n****************", cad);
	unsigned int i=0,j=0,mayor=0;
	//printf("EL SIGUIENTE PASO ES OBTENER LAS PROBABILIDADES\n");
	FILE *fichero2;
	fichero2 = fopen(cad, "a+");
	for(j=0;j<numNodes;j++) {
		if (debugComments) printf("EL NODO %d TIENE %d ENLACES \n",j,ARRAY_DEGREE[j]);
	}
	//el siguiente paso es buscar la mayor cantidad de enlaces por nodo en la red
	for(i=0;i<numNodes;i++) {
		if (ARRAY_DEGREE[i]>mayor) mayor=ARRAY_DEGREE[i];                                          
	}
	//printf("LA MAYOR CANTIDAD DE grado ES %d\n",mayor);  
	j=0;i=0;
	unsigned int cantidad=0;
	for(j=0;j<=mayor;j++) {
		cantidad=0;
		for(i=0;i<numNodes;i++) {
			if (ARRAY_DEGREE[i]==j) cantidad++;
		}
		if(cantidad>0) {
			if (debugComments) printf("CON %d ENLACES HAY %d NODOS\n",j,cantidad);
			fprintf(fichero2, "%d\t%d\n",j,cantidad);
		}
	}
	fprintf(fichero2,"*\n");
	fclose(fichero2);    
	i=0;j=0;
	for(i=0;i<numNodes;i++) {
		if (array_search_1(LISTA_VECINOS[i], ARRAY_DEGREE[i], i)) {
			printf("********************************** Existe un loop en el nodo %d\n", i);
		}                                       
	}
}

/****************************MEAN DEGREE FUNCTION*******************************/
void meanDegreeCalculation (unsigned int connectedNodes, float meanDegree[], unsigned int *ARRAY_DEGREE, unsigned int simulations) {
	unsigned int x=0;
	unsigned int linksSum=0;
	for(x=0;x<connectedNodes;x++) {
		linksSum+=ARRAY_DEGREE[x];		//cantidad de enlaces 
	}
	int ind_prom=connectedNodes/1;
	meanDegree[ind_prom]+=1.0*linksSum/(simulations*connectedNodes);
}
/****************************COPY LINKS FUNCTION*******************************/
int getNumLinks (unsigned int idNode, int AgeNod[], unsigned int connectedNodes, unsigned int *ARRAY_DEGREE, float P, unsigned int initialLinks) {
	float r=random_uniform_f();
	if (debugComments) printf("La probabilidad que el nuevo nodo nazca con kout = 1 es : %f\n", P);
	if (debugComments) printf("El numero aleatorio entre 0 y 1 es: %f\n",r);	
	if (r<P) return initialLinks;
	//if (((AgeNod[idNode]*5+4)<15 || (AgeNod[idNode]*5+4)>50)) return 1;
	//if (r<P) return initialLinks;
	unsigned long random;
	random_uniform(0, connectedNodes-1, &random);
	if (debugComments) printf("El numero aleatorio entre  %d y %d es %ld\n", 0, connectedNodes-1, random);
	if (debugComments) printf("El nodo %d nacio con kout=%d coSIado del nodo %ld\n", idNode, ARRAY_DEGREE[random], random);	
	return ARRAY_DEGREE[random];
}

/****************************PREFERENTIAL LINKING FUNCTIONS*******************************/
void preferentialLinking (unsigned int idNode, int AgeNod[], float P, unsigned int initialLinks,  unsigned int connectedNodes, unsigned int Ai, unsigned int div_Ai, unsigned int *ARRAY_DEGREE, unsigned int *LISTA_VECINOS[]) {
	unsigned int i=0, j=0;
	unsigned long suma_grado_red=0;
	unsigned int numLinks = getNumLinks (idNode, AgeNod, connectedNodes, ARRAY_DEGREE, P, initialLinks);
	//creamos un arreglo auxiliar con los vecinos del nodo idNode
	unsigned int w=0;


	//printf("El nodo %d nacio con %d links\n", idNode, numLinks);
	for(i=0;i<numLinks;i++) {	    
		//el siguiente paso es recorrer los nodos presentes para ver con cual se va a conectar
		unsigned int m=0;
		unsigned long suma_anterior=1,suma=0;
		//el siguiente paso es contar el DEGREE de la red 
		j=0,suma_grado_red=0;
		unsigned int ff=0,acumular=1;
		for(j=0;j<connectedNodes;j++) {
            //if (array_search_1(LISTA_VECINOS[idNode], ARRAY_DEGREE[idNode], j) == 0) {
                if (Ai>0.0) {
                    suma_grado_red+=ARRAY_DEGREE[j]*div_Ai+Ai;
                } else {
                    suma_grado_red+=ARRAY_DEGREE[j];				
                }
            //}
		}
		if (debugComments) printf("El grado de la red es %ld\n", suma_grado_red);
		if (suma_grado_red==0) break; 
		unsigned long random=0;
		random_uniform(1,suma_grado_red,&random);
		if (debugComments) printf("EL NUMERO ALEATORIO entre %d y %ld UNIFORME ES %ld\n",1,suma_grado_red,random);
		for(m=0;m<connectedNodes;m++) {					
            //if (array_search_1(LISTA_VECINOS[idNode], ARRAY_DEGREE[idNode], j) == 0) {
                if (Ai>0.0) {
                    suma+=ARRAY_DEGREE[m]*div_Ai+Ai;
                } else {
                    suma+=ARRAY_DEGREE[m];				
                }				
                if (debugComments) printf("El vector de probabilidad del nodo %d es %ld -> %ld\n",m,suma_anterior,suma);						
                if (random>=suma_anterior && random<=suma) {
                    //if (debugComments) printf("\n\n°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°\n");
                    connectNodes(LISTA_VECINOS, ARRAY_DEGREE, idNode, m);				
                    break;
                }
                suma_anterior=suma+1;
            //}
		}	
	}
}

int getAge ( void ) {
	//Percentage of age in Mexico
	float AgeMex[21]={0.0797321177522273,0.0854220717529027,0.0868438262077878,0.0857578359691141,0.0827058343918928,0.0793007054516409,0.0747601473309034,0.0715815249261463,0.0674812670056469,0.0630280086921119,0.05584721268801,0.0452009849316454,0.0382581386338397,0.0289259630340826,0.0210082966638697,0.0143998417192042,0.00932724757682526,0.00523152089802322,0.00211727228074234,0.000755511148505185,0.0023144886911952};
	float n=0.0, sum=0.0;
	int i=0;
	struct timespec ts1;
	clock_gettime( CLOCK_REALTIME, &ts1 ); 
	n = (rand()%100000000)/100000000.0;
		for (i=0;i<21;i++) {
			sum+=AgeMex[i];
			if (n<sum) return i;
		}
	return i;	
}

void random_uniform(unsigned int M,unsigned long N,unsigned long *numero) {
	*numero = rand () % (N-M+1) + M;//genera un numero aleatorio de M a N
}

float random_uniform_f() {
float n=0.0;
 	struct timespec ts1;
	clock_gettime( CLOCK_REALTIME, &ts1 ); 
	n = (rand()%100000000)/100000000.0;
	//printf("-------> %f\n", n);
return(n);
}


int array_search_1(unsigned int array[], int tamano, int dato) {
	int i=0, senal=0;
	for(i=0;i<tamano;i++) {
		if (dato==array[i]) {
			senal=1;break;
		}
	}
	if(senal==1) {
		return(1);
	} else {
		return(0);
	}
}
