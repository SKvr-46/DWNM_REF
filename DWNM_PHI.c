#include <stdio.h>
#include <complex.h>
#include <math.h>
#include <gsl/gsl_sf_bessel.h>

#define N 20001  // Adjust this to match your data size
#define L M_PI
#define ALPHA 1.0
#define R 1.0
#define Z 1.0
#define N_RANGE 100000

double complex phi_values_omega_domain(double omegas[], int size) {
    double complex phi_values[N];
    
    for (int i = 0; i < size; i++) {
        double omega = omegas[i];
        double complex phi_sum = 0.0 + 0.0 * I;
        
        for (int n = 0; n < N_RANGE; n++) {
            double k = 2 * n * M_PI / L;
            double nu = sqrt((omega / ALPHA) * (omega / ALPHA) - k * k);
            if (nu != nu) {  // check for NaN (imaginary number in sqrt)
                printf("Error: imaginary number in square root\n");
                break;
            }
            if (n == 0) {
                phi_sum += (-I * M_PI / L) * (k / nu) * gsl_sf_bessel_J0(k * R) * cexp(-I * nu * fabs(Z));
            } else {
                phi_sum += (2 * (-I * M_PI / L) * (k / nu) * gsl_sf_bessel_J0(k * R) * cexp(-I * nu * fabs(Z)));
            }
        }
        phi_values[i] = phi_sum;
    }
    return phi_values;
}

void r_Fourier_phi(double complex X[], int size) {  // Fourier inverse transformation
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
    fprintf(gnuplot, "set output 'output_phi.png'\n");
    fprintf(gnuplot, "plot '-' with lines\n");
    
    for(n = 0; n < N; n++) {
        fprintf(gnuplot, "%d %f\n", n, x_r[n]);
    }
    
    fprintf(gnuplot, "e\n");
    fflush(gnuplot);

    pclose(gnuplot);
}

int main() {
    double omegas[N];
    for (int i = 0; i < N; i++) {
        omegas[i] = i - 10000;
    }
    double complex *phi_values = phi_values_omega_domain(omegas, N);
    r_Fourier_phi(phi_values, N);

    return 0;
}
