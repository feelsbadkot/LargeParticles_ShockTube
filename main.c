// ************************************************************* 
// * Численное моделирование процесса течения в ударной трубе. *
// *  Рассматривается осесимметричное течение идеального газа. *
// * Базовая система дифференциальных уравнений: ур-ния Эйлера.*
// *    Используется явная разностная схема метода Давыдова    *
// *           с пересчетом "ЕЕ" по "UE" и "VE"                *
// *************************************************************

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N   550
#define N_D 551
#define M    10
#define M_D  11

double **dmatrix(int nrl, int nrh, int ncl, int nch);
void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch);
void nrerror(char error_text[]);
int main();

int main() {
	// Описание используемых переменных.
  	double R0, A0, RO0, DT, DX, DR, K, 
		   **RO, **U, **V, **P, **E, **RO1,
		   **ROP, **ROB, **UP, **VB,
           **UE, **VE, **EE,
           **RU, **RUU, **RVU, **REU, **RV, **RUV, **RVV, **REV, **RM,
			A1, A2, A3, A4, A5, A6, A7, A8, A10;
  	int I, J, NC;
  	FILE *F01;

    // Определение указателей на переменные.
  	RO = dmatrix(0, N_D, 0, M_D);
  	U = dmatrix(0, N_D, 0, M_D);
  	V = dmatrix(0, N_D, 0, M_D);
  	P = dmatrix(0, N_D, 0, M_D);
  	E = dmatrix(0, N_D, 0, M_D);
  	RO1 = dmatrix(0, N_D, 0, M_D);
  	ROP = dmatrix(0, N_D, 1, M_D);
  	ROB = dmatrix(1, N_D, 0, M_D);
  	UP = dmatrix(0, N_D, 1, M_D);
  	VB = dmatrix(1, N_D, 0, M_D);
  	UE = dmatrix(0, N_D, 0, M_D);
  	VE = dmatrix(0, N_D, 0, M_D);
  	EE = dmatrix(0, N_D, 0, M_D);
  	RU = dmatrix(0, N_D, 0, M_D);
  	RUU = dmatrix(0, N_D, 0, M_D);
  	RVU = dmatrix(0, N_D, 0, M_D);
  	REU = dmatrix(0, N_D, 0, M_D);
  	RV = dmatrix(0, N_D, 0, M_D);
  	RUV = dmatrix(0, N_D, 0, M_D);
  	RVV = dmatrix(0, N_D, 0, M_D);
  	REV = dmatrix(0, N_D, 0, M_D);
  	RM = dmatrix(0, N_D, 0, M_D);

    // Открытие файла для печати расчётной информации.
  	if ((F01 = fopen("./OutMerc02a.txt","w")) == NULL) {
    	printf("File F01 is not open!\n");
    	return 1;
  	}

	// Сообщение о начале счёта.
  	printf("Mercury02, go,go,go!!!\n");
  	fprintf(F01,"Mercury02, go,go,go!!!\n");

    // Константы.
  	R0 = 1; 
	A0 = 340; 
	RO0 = 1.204;
  	DT = 4E-6 * A0 / R0; 
	DX = 0.01 / R0; 
	DR = 0.01 / R0;
  	K = 1.4;

    // Начальные условия.
  	for (I = 1; I <= N; I++) {
		for (J = 1; J <= M; J++) {
	  		if (I <= 20) { 
				RO[I][J] = 20.0 * 1.204 / RO0; 
				P[I][J] = 20.0 * 1.01325E5 / (RO0 * A0 * A0); 
			} else {
				RO[I][J] = 1.204 / RO0; 
				P[I][J] = 1.01325E5 / (RO0 * A0 * A0); 
			}
	  		U[I][J] = 0.0; 
			V[I][J] = 0.0;
	  		E[I][J] = P[I][J] / (RO[I][J] * (K - 1));
		}
  	}

    //***************************************
    //*          Цикл по времени.           *
    //***************************************

  	for (NC = 1; NC < 100000; NC++) {
		// ГУ: правая граница.
		for (J = 1; J <= M; J++) {
			RO[N+1][J] = RO[N][J]; 
			U[N+1][J] = -U[N][J]; 
			V[N+1][J] = V[N][J]; 
			P[N+1][J] = P[N][J];
		}

    	// ГУ: нижняя граница - ось симметрии.
		for (I = 1; I <= N; I++) {
			RO[I][0] = RO[I][1]; 
			U[I][0] = U[I][1]; 
			V[I][0] = -V[I][1]; 
			P[I][0] = P[I][1]; 
		}

        // ГУ: левая граница.
		for (J = 1; J <= M; J++) {
			RO[0][J] = RO[1][J]; 
			U[0][J] = -U[1][J]; 
			V[0][J] = V[1][J];
			P[0][J] = P[1][J]; 
		}

	    // ГУ: верхняя граница.
		for (I = 1; I <= N; I++) {
			RO[I][M+1] = RO[I][M]; 
			U[I][M+1] = U[I][M]; 
			V[I][M+1] = -V[I][M]; 
			P[I][M+1] = P[I][M];
		}

	    //********************************
    	//"Эйлеров" этап метода Давыдова.*
	    //********************************

	    // Вычисление "эйлеровых" скоростей.
    	for (I = 1; I <= N; I++) {
      		for (J = 1; J <= M; J++) {
        		A1 = (P[I][J] + P[I+1][J]) / 2; 
				A2 = (P[I-1][J] + P[I][J]) / 2;
        		A3 = (P[I][J] + P[I][J+1]) / 2; 
				A4 = (P[I][J-1] + P[I][J]) / 2;
        		UE[I][J] = U[I][J] - (A1 - A2) * DT / (RO[I][J] * DX);
        		VE[I][J] = V[I][J] - (A3 - A4) *DT / (RO[I][J] * DR);
      		}
    	}

        // Е.ГУ: правая граница.
    	for (J = 1; J <= M; J++) { 
			UE[N+1][J] = -UE[N][J]; 
			VE[N+1][J] = VE[N][J]; 
		}

        // Е.ГУ: нижняя граница - ось симметрии.
    	for (I = 1; I <= N; I++) { 
			UE[I][0] = UE[I][1]; 
			VE[I][0] = -VE[I][1]; 
		}

        // Е.ГУ: левая граница.
    	for (J = 1; J <= M; J++) { 
			UE[0][J] = -UE[1][J];
			VE[0][J] = VE[1][J]; 
		}

	    // Е.ГУ: верхняя граница.
		for (I = 1; I <= N; I++) { 
			UE[I][M+1] = UE[I][M]; 
			VE[I][M+1] = -VE[I][M]; 
		}

	    // Вычисление "эйлеровой" полной удельной энергии.	
		for (I = 1; I <= N; I++) {
			for (J = 1; J <= M; J++) {
				A1 = (P[I][J] + P[I+1][J]) / 2; 
				A2 = (P[I-1][J] + P[I][J]) / 2;
				A3 = (P[I][J] + P[I][J+1]) / 2; 
				A4 = (P[I][J-1] + P[I][J]) / 2;
				A5 = (UE[I][J] + UE[I+1][J]) / 2; 
				A6 = (UE[I-1][J] + UE[I][J]) / 2;
				A7 = (VE[I][J] + VE[I][J+1]) / 2; 
				A8 = (VE[I][J-1] + VE[I][J]) / 2;
				EE[I][J] = E[I][J]- (A1 * A5 - A2 * A6) * DT / (RO[I][J] * DX) -
						   (J * A3 * A7 - (J - 1) * A4 * A8) * DT / (RO[I][J] * (J - 0.5) * DR);
			}
		}

        // ЕE.ГУ: правая граница.
    	for (J = 1; J <= M; J++) { EE[N+1][J] = EE[N][J]; }

        // ЕE.ГУ: нижняя граница - ось симметрии.
    	for (I = 1; I <= N; I++) { EE[I][0] = EE[I][1]; }

        // ЕЕ.ГУ: левая граница.
    	for (J = 1; J <= M; J++) { EE[0][J] = EE[1][J]; }

	    // ЕЕ.ГУ: верхняя граница.
    	for (I = 1; I <= N; I++) { EE[I][M+1] = EE[I][M]; }

	    // Анализ устойчивости вычислений.
		for (I = 1; I <= N; I++) {
			for (J = 1; J <= M; J++) {
				if (P[I][J] < 0) {
					printf("Attention!!! P[I][J]= %e %d %d %d\n", P[I][J], I, J, NC);
					fprintf(F01, "Attention!!! P[I][J]= %e %d %d %d\n", P[I][J], I, J, NC);
					fprintf(F01, "     Parameter PR[I][J]=\n");

					for (I = 0; I <= N + 1; I++) {
						fprintf(F01, "%3d ", I);
						for (J=0; J <= M+1; J++) { 
							RM[I][J]=P[I][J] * A0 * A0 * RO0 / 1E6;
							fprintf(F01, "%8.3f", RM[I][J]);
						}
						fprintf(F01,"\n");
					}

					fprintf(F01,"     Parameter U1R[I][J]=\n");
					for (I = 0; I <= N + 1; I++) {
						fprintf(F01, "%3d ", I);
						for (J = 0; J <= M+1; J++) { 
							RM[I][J] = U[I][J] * A0;
							fprintf(F01, "%8.2f", RM[I][J]);
						}
						fprintf(F01, "\n");
					
					}
					return 20;
				}
			}
		}

		//**********************************
		//"Лагранжев" этап метода Давыдова.*
		//**********************************
		for (I = 0; I <= N; I++) {
			for (J = 0; J <= M; J++) {
				// Вдоль оси "0Х".
				A1 = (UE[I][J] + UE[I+1][J]) / 2;
				if (A1 >= 0) {
					RU[I][J] = RO[I][J] * A1;
					RUU[I][J] = RU[I][J] * UE[I][J];
					RVU[I][J] = RU[I][J] * VE[I][J];
					REU[I][J] = RU[I][J] * EE[I][J];
				} else {
					RU[I][J] = RO[I+1][J] * A1;
					RUU[I][J] = RU[I][J] * UE[I+1][J];
					RVU[I][J] = RU[I][J] * VE[I+1][J];
					REU[I][J] = RU[I][J] * EE[I+1][J];
				}

				// Вдоль оси "0R".
				A1 = (VE[I][J] + VE[I][J+1]) / 2;
				if (A1 >= 0) {
					RV[I][J] = RO[I][J] * A1;
					RUV[I][J] = RV[I][J] * UE[I][J];
					RVV[I][J] = RV[I][J] * VE[I][J];
					REV[I][J] = RV[I][J] * EE[I][J];
				} else {
					RV[I][J] = RO[I][J+1] * A1;
					RUV[I][J] = RV[I][J] * UE[I][J+1];
					RVV[I][J] = RV[I][J] * VE[I][J+1];
					REV[I][J] = RV[I][J] * EE[I][J+1];
				}
			}
		}

		//***************************************
		//"Заключительный" этап метода Давыдова.*
		//***************************************
		for (I = 1; I <= N; I++) {
			for (J = 1; J <= M; J++) {
				RO1[I][J] = RO[I][J] - (RU[I][J] - RU[I-1][J]) * DT / DX-
							(J * RV[I][J] - (J - 1) * RV[I][J-1]) * DT / ((J - 0.5) * DR);
				A1 = DT / RO1[I][J];
				U[I][J] = RO[I][J] * UE[I][J] / RO1[I][J]-
						  (RUU[I][J] - RUU[I-1][J]) * A1 / DX -
						  (J * RUV[I][J] - (J - 1) * RUV[I][J-1]) * A1 / ((J - 0.5) * DR);
				V[I][J] = RO[I][J] * VE[I][J] / RO1[I][J]-
						  (RVU[I][J] - RVU[I-1][J]) * A1 / DX -
						  (J * RVV[I][J] - (J - 1) * RVV[I][J-1]) * A1 / ((J - 0.5) * DR);
				E[I][J] = RO[I][J] * EE[I][J] / RO1[I][J]-
						  (REU[I][J] - REU[I-1][J]) * A1 / DX - 
						  (J * REV[I][J] - (J - 1) * REV[I][J-1]) * A1 / ((J - 0.5) * DR);
			}
		}

		// Дополнительные параметры
		for (I = 1; I <= N; I++) {
			for (J=1; J<=M; J++) { 
				P[I][J] = (K - 1) * RO1[I][J] * (E[I][J] - (U[I][J] * U[I][J] + V[I][J] * V[I][J]) / 2);
			}
		}
		
		// Переприсвоение значений.
		for (I = 1; I <= N; I++) {
			for (J = 1; J <= M; J++) {
				RO[I][J] = RO1[I][J];
			}
		}

		// Печать результатов расчёта.

		// Печать текущая.
		if (!(NC % 25)) {
			A10 = A0 * A0 * RO0 / 1E6;
			A1 = P[1][5] * A10; 
			A2 = P[20][5] * A10;
			A3 = P[40][5] * A10;
			A4 = P[60][5] * A10; 
			A5 = P[80][5] * A10; 
			A6 = U[100][5] * A0;
			fprintf(F01,"%5d %8.3f %8.3f %8.3f %8.3f %8.3f %8.2f\n",
					NC, A1, A2, A3, A4, A5, A6);
		}

		// Печать параметров по полному полю течения.
		if (!(NC % 100)) {
			fprintf(F01, "Parameters of stream %d\n", NC);
			fprintf(F01, "     Parameter PR[I][J]=\n");
			for (I = 1; I <= N; I++) {
				fprintf(F01, "%3d ", I);
				for (J = 1; J <= M; J++) { 
					RM[I][J] = P[I][J] * A0 * A0 * RO0 / 1E6;
					fprintf(F01, "%8.3f", RM[I][J]);
				}
				fprintf(F01, "\n");
			}
			fprintf(F01, "     Parameter UR[I][J]=\n");
			for (I = 1; I <= N; I++) {
				fprintf(F01, "%3d ", I);
				for (J = 1; J <= M; J++) { 
					RM[I][J] = U[I][J] * A0;
					fprintf(F01, "%8.2f", RM[I][J]);
				}
				fprintf(F01, "\n");
			}
		} 
		if (NC == 1000)  {
			fclose(F01);
			printf("The END. Good time!\n");
			return 1000;
		}
  	}  // Конец цикла по времени!

    // Освобождение динамической памяти для переменных.
  	free_dmatrix(RO, 0, N_D, 0, M_D);
  	free_dmatrix(U, 0, N_D, 0, M_D);
  	free_dmatrix(V, 0, N_D, 0, M_D);
  	free_dmatrix(P, 0, N_D, 0, M_D);
  	free_dmatrix(E, 0, N_D, 0, M_D);
  	free_dmatrix(RO1, 0, N_D, 0, M_D);
  	free_dmatrix(ROP, 0, N_D, 1, M_D);
  	free_dmatrix(ROB, 1, N_D, 0, M_D);
  	free_dmatrix(UP, 0, N_D, 1, M_D);
  	free_dmatrix(VB, 1, N_D, 0, M_D);
  	free_dmatrix(UE, 0, N_D, 0, M_D);
  	free_dmatrix(VE, 0, N_D, 0, M_D);
  	free_dmatrix(EE, 0, N_D, 0, M_D);
  	free_dmatrix(RU, 0, N_D, 0, M_D);
  	free_dmatrix(RUU, 0, N_D, 0, M_D);
  	free_dmatrix(RVU, 0, N_D, 0, M_D);
  	free_dmatrix(REU, 0, N_D, 0, M_D);
  	free_dmatrix(RV, 0, N_D, 0, M_D);
  	free_dmatrix(RUV, 0, N_D, 0, M_D);
  	free_dmatrix(RVV, 0, N_D, 0, M_D);
  	free_dmatrix(REV, 0, N_D, 0, M_D);

  	return 0;
} // Конец Mercury02 !!!

// Двумерная матрица удвоенной точности
double **dmatrix(int nrl, int nrh, int ncl, int nch) {
	int i;
	double **m;
	m = (double**)calloc((unsigned long) (nrh - nrl + 1), sizeof(double*));
	if (!m) { nrerror("allocation failure 1 in dmatrix()"); }
	m -= nrl;
	for (i = nrl; i <= nrh; i++) {
		m[i] = (double*)calloc((unsigned long)(nch - ncl + 1), sizeof(double));
		if (!m[i]) { nrerror("allocation failure 2 in dmatrix()"); }
		m[i] -= ncl;
	}
	return m;
}

void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch) {
	int i;
	for(i = nrh; i >= nrl; i--) { free((char*)(m[i] + ncl)); }
	free((char*)(m + nrl));
}

void nrerror(char error_text[]) {
	void exit(); 
	fprintf(stderr, "Numerical Recipes run-time error...\n");
	fprintf(stderr, "%s\n", error_text);
	fprintf(stderr, "...now exiting to system...\n");
	exit(1);
}
