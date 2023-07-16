#include <stdio.h>
#include <complex.h>
#include <math.h>

#define N 20001  // Adjust this to match your data size
#define L 25.0
#define ALPHA 5.0
#define R 10.0
#define Z 1.0
#define N_RANGE 1000
#define f_s 500.0
#define t_fin 20.0
#define dt 1/f_s
#define OMEGA_IMAG (-1j * 2 * M_PI / 1000)

double t[N];
double omegas[N];

double complex* u_values_omega_domain(double complex omegas[], int size) {
    double complex* u_values = malloc(size * sizeof(double complex));
    
    for (int i = 0; i < size; i++) {
        double complex omega = omegas[i];
        double complex u_sum = 0.0 + 0.0 * I;
        
        for (int n = 0; n < N_RANGE; n++) {
            double k = 2 * n * M_PI / L;
            double complex nu_val = (omega / ALPHA) * (omega / ALPHA) - k * k;
            double complex nu = cimag(sqrt(nu_val)) < 0 ? sqrt(nu_val) : -1 * sqrt(nu_val);
            if (cimag(nu) > 0) {
                printf("Error: imaginary number in square root\n");
                break;
            }
            if (n == 0) {
                u_sum += (I * M_PI / L) * (k * k / nu) * j1(k * R) * cexp(-I * nu * fabs(Z));
            } else {
                u_sum += (2 * (I * M_PI / L) * (k * k / nu) * j1(k * R) * cexp(-I * nu * fabs(Z)));
            }
        }
        u_values[i] = u_sum;
    }
    return u_values;
}

void r_Fourier(double complex X[], int size) {  // Fourier inverse transformation
    int n, k;
    double x_r[N];
    
    for(n = 0; n < size; n++) {
        double complex r = 0.0 + 0.0 * I;
            
        for(k = 0; k < size; k++) {
            r += X[k] * cexp((2 * M_PI * k * n * I) / size);
        }
            
        x_r[n] = creal(r) / size;
    }
    
    FILE *gnuplot = popen("gnuplot", "w");
    fprintf(gnuplot, "set term png\n");
    fprintf(gnuplot, "set output 'output.png'\n");
    fprintf(gnuplot, "set xlabel 'time t'\n");
    fprintf(gnuplot, "set ylabel 'u'\n");
    fprintf(gnuplot, "plot '-' with lines\n");
    
    for(n = 0; n < N; n++) {
        fprintf(gnuplot, "%f %f\n", t[n], x_r[n]);
    }
    
    fprintf(gnuplot, "e\n");
    fflush(gnuplot);

    pclose(gnuplot);
}

int main() {
    // Initialize the time array and frequency array
    for (int i = 0; i < N; i++) {
        t[i] = i * dt;
        omegas[i] = (i/N)*f_s - f_s/2;
    }

    double complex* u_values = u_values_omega_domain(omegas, N);
    r_Fourier(u_values, N);

    // Don't forget to free the dynamically allocated memory
    free(u_values);

    return 0;
}
