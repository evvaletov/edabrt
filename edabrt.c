/****************************************************************************
 *              edabrt: Electrostatic Deflector Aberrations                 *
 *                        E. Valetov & M. Berz                              *
 *                         Created 25-Jan-2018                              *
 *                       Email: valetove@msu.edu                            *
 *                                                                          *
 * This program computes the first and second order aberrations of an       *
 * electrostatic deflector in the x-a plane using exact analytic formulas.  *
 * Please refer to the following report for a derivation of these formulas: *
 * E. Valetov and M. Berz                                                   *
 * Derivation of Analytic Formulas for Electrostatic Deflector Aberrations, *
 * and Comparison with the Code COSY INFINITY                               *
 * MSUHEP-180212, Michigan State University (2018)                          *
 ****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

double calculate_aberration(double radius, double ang, double n1, double n2,
                            int index1, int index2) {
    double s = radius * (M_PI / 180) * ang; // arc length
    double h = 1 / radius; // curvature
    if (n1 < 3) {
        double s3n1 = h * sqrt(3 - n1);
        switch (index1) {
            case 1:
                switch (index2) {
                    case 10:
                        return cos(s3n1 * s);
                    case 01:
                        return sin(s3n1 * s) / s3n1;
                    case 20:
                        return -4 * h * (9 * n1 + 2 * n2 - 15 +
                                         (6 * n1 + n2 - 12) * cos(s3n1 * s)) *
                               pow(sin(s3n1 * s / 2), 2) / (3 * (n1 - 3));
                    case 11:
                        return -2 * pow(h, 3) * (3 - 3 * n1 - n2 +
                                                 (6 * n1 + n2 - 12) *
                                                 cos(s3n1 * s)) *
                               sin(s3n1 * s) /
                               (3 * pow(h * h * (3 - n1), 1.5));
                    case 02:
                        return -4 * (3 - 3 * n1 - n2 +
                                     (6 * n1 + n2 - 12) * cos(s3n1 * s)) *
                               pow(sin(s3n1 * s / 2), 2) /
                               (3 * h * pow((n1 - 3), 2));
                    default:
                        return 0;
                }
            case 2:
                switch (index2) {
                    case 10:
                        return -s3n1 * sin(s3n1 * s);
                    case 01:
                        return cos(s3n1 * s);
                    case 20:
                        return 2 * pow(h, 3) * (3 * n1 + n2 - 3) *
                               (sin(s3n1 * s) + sin(2 * s3n1 * s)) /
                               (3 * s3n1);
                    case 11:
                        return -4 * h * (3 * n1 + n2 - 3) *
                               (1 + 2 * cos(s3n1 * s)) *
                               pow(sin(s3n1 * s / 2), 2) /
                               (3 * (n1 - 3));
                    case 02:
                        return -2 * pow(h, 3) * (15 - 9 * n1 - 2 * n2 +
                                                 2 * (3 * n1 + n2 - 3) *
                                                 cos(s3n1 * s)) *
                               sin(s3n1 * s) /
                               (3 * pow(h * h * (3 - n1), 1.5));
                    default:
                        return 0;
                }
            default:
                return 0;
        }
    } else if (n1 > 3) {
        double s3n1 = h * sqrt(n1 - 3);
        switch (index1) {
            case 1:
                switch (index2) {
                    case 10:
                        return cosh(s3n1 * s);
                    case 01:
                        return -sinh(s3n1 * s) / s3n1;
                    case 20:
                        return h * (4 * (3 * n1 + n2 - 3) * cosh(s3n1 * s) -
                                    3 * (4 * n1 + n2 - 6) *
                                    (2 * cosh(2 * s3n1 * s) - 1) -
                                    (6 + n2) * cosh(4 * s3n1 * s)) /
                               (6 * (n1 - 3));
                    case 11:
                        return (4 * (3 * n1 + n2 - 3) * sinh(s3n1 * s) -
                                (6 + n2) * sinh(4 * s3n1 * s)) /
                               (6 * pow(n1 - 3, 1.5));
                    case 02:
                        return -(4 * (9 * n1 + 2 * n2 - 15) * cosh(s3n1 * s) -
                                 3 * (4 * n1 + n2 - 6) *
                                 (2 * cosh(2 * s3n1 * s) + 1) +
                                 (6 + n2) * cosh(4 * s3n1 * s)) /
                               (6 * h * pow(n1 - 3, 2));
                    default:
                        return 0;
                }
            case 2:
                switch (index2) {
                    case 10:
                        return s3n1 * sinh(s3n1 * s);
                    case 01:
                        return cosh(s3n1 * s);
                    case 20:
                        return 2 * h * h * h * (3 * n1 + n2 - 3) *
                               (sinh(s3n1 * s) + sinh(2 * s3n1 * s)) /
                               (3 * s3n1);
                    case 11:
                        return 2 * h * (3 * n1 + n2 - 3) *
                               (cosh(2 * s3n1 * s) - cosh(s3n1 * s)) /
                               (3 * (n1 - 3));
                    case 02:
                        return 2 * (15 - 9 * n1 - 2 * n2 +
                                    2 * (3 * n1 + n2 - 3) *
                                    cosh(s3n1 * s)) *
                               sinh(s3n1 * s) /
                               (3 * pow(n1 - 3, 1.5));
                    default:
                        return 0;
                }
            default:
                return 0;
        }
    } else {
        switch (index1) {
            case 1:
                switch (index2) {
                    case 10:
                        return 1;
                    case 01:
                        return s;
                    case 20:
                        return h * h * h * (6 + n2) * s * s;
                    case 11:
                        return 2 * h * s +
                               h * h * h * (6 + n2) * s * s * s / 3;
                    case 02:
                        return h * s * s * (6 + h * h * (6 + n2) * s * s) / 6;
                    default:
                        return 0;
                }
            case 2:
                switch (index2) {
                    case 10:
                        return 0;
                    case 01:
                        return 1;
                    case 20:
                        return 2 * h * h * h * (6 + n2) * s;
                    case 11:
                        return h * h * h * (6 + n2) * s * s;
                    case 02:
                        return 2 * h * s * (h * h * (6 + n2) * s * s - 3) / 3;
                    default:
                        return 0;
                }
            default:
                return 0;
        }
    }
}

int print_aberrations(double radius, double ang, double n1, double n2,
                      int index1) {
    double aberration;
    int max_order = 2;
    int number_of_variables = 2;
    int first_non_zero = 0;
    char index_string[20];
    int counter = 1;
    for (int order = 1; order <= max_order; order++) {
        int index_array[number_of_variables + 1];
        for (int j = 0; j <= number_of_variables; j++)
            index_array[j] = 0;
        index_array[0] = order;
        while (index_array[number_of_variables] == 0) {
            int index2 = 0;
            for (int j = 0; j <= number_of_variables - 1; j++)
                index2 = index2 * 10 + index_array[j];
            aberration = calculate_aberration(radius, ang, n1, n2, index1,
                                              index2);
            sprintf(index_string, "%i", index_array[0]);
            if (number_of_variables > 1)
                for (int j = 1; j <= number_of_variables - 1; j++) {
                    snprintf(index_string, sizeof index_string, "%s %i",
                             index_string, index_array[j]);
                }
            if (aberration != 0) {
                if (counter == 1)
                    printf("     I  COEFFICIENT           ORDER EXPONENTS\n");
                printf("     %i % 16.15le   %i   %s\n", counter, aberration,
                       order, index_string);
                counter++;
            }
            for (int j = 0; j <= number_of_variables - 1; j++) {
                if (index_array[j] != 0) {
                    first_non_zero = j;
                    break;
                }
            }
            if (index_array[first_non_zero] == order) {
                index_array[0] = order - 1;
                index_array[first_non_zero + 1] = 1;
            } else {
                index_array[first_non_zero]--;
                index_array[first_non_zero + 1]++;
            }
        }
    }
    if (counter == 1)
        printf("     ALL COMPONENTS ZERO\n");
    printf("     --------------------------------------\n");
}

int invalid_option_exit(char *option) {
    printf("edabrt: invalid option -- %s\n", option);
    printf("Try 'edabrt --help' for more information.\n");
    exit(-1);
}

int main(int argc, char **argv) {
    double r; // reference radius
    double ang; // central angle
    double n1, n2; // first and second order inhomogeneity coefficients
    double dummy;
    int read_code;
    printf("----------------------------------------------------------\n");
    printf("      edabrt: Electrostatic Deflector Aberrations         \n");
    printf("                 E. Valetov & M. Berz                     \n");
    printf("                  Created 25-Jan-2018                     \n");
    printf("                Email: valetove@msu.edu                   \n");
    printf("----------------------------------------------------------\n");
    if (argc >= 2)
        if (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "/h") == 0 ||
            strcmp(argv[1], "--help") == 0) {
            printf("\nThis program computers the first and second order ");
            printf("aberrations of an electrostatic deflector\n");
            printf("in the horizontal x-a plane using exact analytic ");
            printf("formulas.\n");
            printf("\nINTERACTIVE MODE\n");
            printf("Run the program and follow the prompts to specify the ");
            printf("electrostatic deflector parameters.\n");
            printf("\nCOMMAND-LINE ARGUMENTS\n");
            printf("Electrostatic deflector parameters may be optionally ");
            printf("supplied using the command line:\n");
            printf("edabrt [r ang n1 n2] [--help]\n");
            printf("    r       Reference orbit radius in meters\n");
            printf("    ang     Central angle spanning the deflector ");
            printf("in degrees\n");
            printf("    n1      First order electrostatic field ");
            printf("inhomogeneity coefficient\n");
            printf("    n2      Second order electrostatic field ");
            printf("inhomogeneity coefficient\n");
            printf("    --help  This information\n");
            exit(0);
        }
    switch (argc) {
        case 5:
            if (sscanf(argv[1], "%lf", &r) != 1) invalid_option_exit(argv[1]);
            if (r <= 0) {
                printf("\nedabrt: supplied radius r is not positive\n");
                exit(-1);
            }
            if (sscanf(argv[2], "%lf", &ang) != 1)
                invalid_option_exit(argv[2]);
            if (sscanf(argv[3], "%lf", &n1) != 1)
                invalid_option_exit(argv[3]);
            if (sscanf(argv[4], "%lf", &n2) != 1)
                invalid_option_exit(argv[4]);
            printf("\nReference radius r = % 16.15le m\n", r);
            printf("Central angle ang = % 16.15le°\n", ang);
            printf("1st order inhomogeneity coefficient n1 = % 16.15le\n",
                   n1);
            printf("2nd order inhomogeneity coefficient n2 = % 16.15le\n",
                   n2);
            break;
        case 1:
            printf("\n");
            do {
                fseek(stdin, 0, SEEK_END);
                printf("Please enter the reference orbit radius r ");
                printf("in [m].\n> ");
                stdin = freopen(NULL, "r", stdin);
                read_code = scanf("%lf", &r);
                if (read_code != 1)
                    printf("Not a numerical value.\n");
                else if (r <= 0) {
                    printf("The radius must be positive.\n");
                    read_code = 0;
                }
            } while (read_code != 1);
            do {
                printf("Please enter the central angle ang spanning the ");
                printf("deflector in [°].\n> ");
                stdin = freopen(NULL, "r", stdin);
                read_code = scanf("%lf", &ang);
                if (read_code != 1)
                    printf("Not a numerical value.\n");
            } while (read_code != 1);
            do {
                printf("Please enter the first order inhomogeneity ");
                printf("coefficient n1.\n> ");
                stdin = freopen(NULL, "r", stdin);
                read_code = scanf("%lf", &n1);
                if (read_code != 1)
                    printf("Not a numerical value.\n");
            } while (read_code != 1);
            do {
                printf("Please enter the first order inhomogeneity ");
                printf("coefficient n2.\n> ");
                stdin = freopen(NULL, "r", stdin);
                read_code = scanf("%lf", &n2);
                if (read_code != 1)
                    printf("Not a numerical value.\n");
            } while (read_code != 1);
            break;
        default:
            printf("\n");
            for (int j = 1; j < argc; j++) {
                if (sscanf(argv[j], "%lf", &dummy) != 1) {
                    printf("edabrt: invalid option -- %s\n", argv[j]);
                    printf("Try 'edabrt --help' for more information.\n");
                    exit(-1);
                }
            }
            printf("edabrt: 4 numerical arguments expected, %i supplied\n",
                   argc - 1);
            printf("Try 'edabrt --help' for more information.\n");
            exit(-1);
    }
    printf("\nFirst and second order aberrations in the x-a plane:\n\n");
    printf("(x|...)\n");
    print_aberrations(r, ang, n1, n2, 1);
    printf("(a|...)\n");
    print_aberrations(r, ang, n1, n2, 2);
}