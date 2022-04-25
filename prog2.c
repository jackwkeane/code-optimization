/*
    First and Last name: Jack Keane
    Lehigh email address: jwk324@lehigh.edu
*/

#include <stdio.h>
#include <math.h>
#include <time.h> /* Included for the use of clock */
#include <stdlib.h>
#define PI 3.14159265

#define IDENT 0 /* Constant value used for loop unrolling */
#define OP + /* Constant operation used for loop unrolling */ 

double realOut_ref[N][N] = {0}, imagOut_ref[N][N] = {0}, amplitudeOut_ref[N][N] = {0}; /* Global arrays for reference testing */

void fft_orig(); /* fft function - original */
void fft_test(); /* fft function - used for checking testing correctness */
void fft_code_motion(); /* fft function - optimization: code motion */
void fft_elim_mem_refs(); /* fft function - optimization: eliminated memory references when possible */
void fft_loop_unrolling_2_1(); /* fft functon - optimization: 2x1 loop unrolling */
void fft_loop_unrolling_2_2(); /* fft function - optimization: 2x2 loop unrolling */
void fft_loop_unrolling_4_1(); /* fft function - optimization: 4x1 loop unrolling */
void fft_loop_unrolling_8_1(); /* fft function - optimization: 8x1 loop unrolling */

void check_results(double[N][N], double[N][N], double[N][N]); /* check results - function used when TEST is defined. */
int floating_point_equal(double a, double b); /* check to see if floating points are equal - function used when TEST is defined */

int main() { 
    #ifndef TEST
    /* Orignal FFT */
    printf("%-15ld\t", N);
    clock_t time = clock();
    fft_orig();
    time = clock() - time;
    printf("%-15ld\t", time);
    /* FFT w/ code motion */
    time = 0;
    time = clock();
    fft_code_motion();
    time = clock() - time;
    printf("%-15ld\t", time);
    /* FFT w/ memory references elminated when possible */
    time = 0;
    time = clock();
    fft_elim_mem_refs();
    time = clock() - time;
    printf("%-15ld\t", time);
    /* FFT w/ 2x1 loop unrolling */
    time = 0;
    time = clock();
    fft_loop_unrolling_2_1();
    time = clock() - time;
    printf("%-15ld\t", time);
    /* FFT w/ 2x2 loop unrolling */
    time = 0;
    time = clock();
    fft_loop_unrolling_2_2();
    time = clock() - time;
    printf("%-15ld\t", time);
    /* FFT w/ 4x1 loop unrolling */
    time = 0;
    time = clock();
    fft_loop_unrolling_4_1();
    time = clock() - time;
    printf("%-15ld\t", time);
    /* FFT w/ 8x1 loop unrolling */
    time = 0;
    time = clock();
    fft_loop_unrolling_8_1();
    time = clock() - time;
    printf("%-15ld\n", time);
    #else
    /* Print formatting */
    printf("---------Testing for N = %d---------\nCode\t\tNo-mem\t\t2x1\t\t2x2\t\t4x1\t\t8x1\n", N);
    /* Performing functions w TEST defined */
    fft_orig(); 
    fft_code_motion();
    fft_elim_mem_refs();
    fft_loop_unrolling_2_1();
    fft_loop_unrolling_2_2();
    fft_loop_unrolling_4_1();
    fft_loop_unrolling_8_1();
    /* Print formatting */
    printf("\n------------------------------------\n\n");
    #endif
    return 0;
}

/* FFT function - original function */
void fft_orig(){
    double realOut[N][N] = {0}, imagOut[N][N] = {0}, inputData[N][N] = {0}, amplitudeOut[N][N] = {0};
    int height = N,  width = N;
    int yWave, xWave, ySpace, xSpace, i, j;
   //initialize the array inputData
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
            #ifndef TEST
            inputData[i][j] = rand();
            #else
            inputData[i][j] = 1.0;
            #endif
    // Two outer loops iterate on output data.
    for (yWave = 0; yWave < height; yWave++) {
        for (xWave = 0; xWave < width; xWave++) {
            // Two inner loops iterate on input data.
            for (ySpace = 0; ySpace < height; ySpace++) {
                for (xSpace = 0; xSpace < width; xSpace++) {
                    // Compute real, imag, and ampltude.
                    realOut[yWave][xWave] += (inputData[ySpace][xSpace] * cos(2 * PI * ((1.0 * xWave * xSpace / width) + (1.0 * yWave * ySpace / height)))) / sqrt(width * height);
                    imagOut[yWave][xWave] -= (inputData[ySpace][xSpace] * sin(2 * PI * ((1.0 * xWave * xSpace / width) + (1.0 * yWave * ySpace / height)))) / sqrt(width * height);
                    amplitudeOut[yWave][xWave] = sqrt(realOut[yWave][xWave] * realOut[yWave][xWave]+ imagOut[yWave][xWave] * imagOut[yWave][xWave]);
                }
            }
        }
    }

    #ifdef TEST
    for (int x = 0; x < N; x++) {
        for (int y = 0; y < N; y++) {
            realOut_ref[x][y] = realOut[x][y];
            imagOut_ref[x][y] = imagOut[x][y];
            amplitudeOut_ref[x][y] = amplitudeOut[x][y];
        }
    }
    #endif


}

/* FFT test function - duplicate function of fft_orig that was used for testing the correctness of the check_results() function. */ 
void fft_test() {
    double realOut[N][N] = {0}, imagOut[N][N] = {0}, inputData[N][N] = {0}, amplitudeOut[N][N] = {0};
    int height = N,  width = N;
    int yWave, xWave, ySpace, xSpace, i, j;
   //initialize the array inputData
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
            #ifndef TEST
            inputData[i][j] = rand();
            #else
            inputData[i][j] = 1.0;
            #endif
    // Two outer loops iterate on output data.
    for (yWave = 0; yWave < height; yWave++) {
        for (xWave = 0; xWave < width; xWave++) {
            // Two inner loops iterate on input data.
            for (ySpace = 0; ySpace < height; ySpace++) {
                for (xSpace = 0; xSpace < width; xSpace++) {
                    // Compute real, imag, and ampltude.
                    realOut[yWave][xWave] += (inputData[ySpace][xSpace] * cos(2 * PI * ((1.0 * xWave * xSpace / width) + (1.0 * yWave * ySpace / height)))) / sqrt(width * height);
                    imagOut[yWave][xWave] -= (inputData[ySpace][xSpace] * sin(2 * PI * ((1.0 * xWave * xSpace / width) + (1.0 * yWave * ySpace / height)))) / sqrt(width * height);
                    amplitudeOut[yWave][xWave] = sqrt(realOut[yWave][xWave] * realOut[yWave][xWave]+ imagOut[yWave][xWave] * imagOut[yWave][xWave]);
                }
            }
        }
    }
    
    check_results(realOut, imagOut, amplitudeOut);
}

/* FFT function – optimized using code motion. 
    Summary of changes: 
        Created variables, demon_cm, yWave_cm, xWave_cm, and two_pi to reduce the total amount of calculations neccessary.
        Moved the assignment for each index of amplitudeOut outside of the inner most for-loop and placed it after the 3rd for-loop. 

*/
void fft_code_motion() {
    double realOut[N][N] = {0}, imagOut[N][N] = {0}, inputData[N][N] = {0}, amplitudeOut[N][N] = {0};
    int height = N,  width = N;
    int yWave, xWave, ySpace, xSpace, i, j;
   //initialize the array inputData
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
            #ifndef TEST
            inputData[i][j] = rand();
            #else
            inputData[i][j] = 1.0;
            #endif

    /* Calculation moved outside the inner most loop - using denomintar variable to reduce the amount of times sqrt iscalled. */
    double denom_cm = sqrt(width * height); 
    double two_pi = 2 * PI; /* Saving the program from having to repeat calculations */
            

    // Two outer loops iterate on output data.
    for (yWave = 0; yWave < height; yWave++) {
        /* Calculations replaced from the inner three loops */
        double yWave_cm = 1.0 * yWave / height; 

        for (xWave = 0; xWave < width; xWave++) {
            // Two inner loops iterate on input data.
        
            /* Calculations replaced from the inner two loops */
            double xWave_cm = 1.0 * xWave / width; 

            for (ySpace = 0; ySpace < height; ySpace++) {
                for (xSpace = 0; xSpace < width; xSpace++) {
                    // Compute real, imag, and ampltude.
                    realOut[yWave][xWave] += (inputData[ySpace][xSpace] * cos(two_pi * ((1.0 * xWave_cm * xSpace) + (1.0 * yWave_cm * ySpace)))) / denom_cm;
                    imagOut[yWave][xWave] -= (inputData[ySpace][xSpace] * sin(two_pi * ((1.0 * xWave_cm * xSpace) + (1.0 * yWave_cm * ySpace)))) / denom_cm;
                }
            }

            amplitudeOut[yWave][xWave] = sqrt(realOut[yWave][xWave] * realOut[yWave][xWave]+ imagOut[yWave][xWave]* imagOut[yWave][xWave]);  /* Relocated outside and after the loop */
        }
    }

    /* Checking to see if code matches */
    #ifdef TEST
    check_results(realOut, imagOut, amplitudeOut);
    #endif
}

/* FFT function – eliminated memory references when possible. 
    Summary of changes:
        Prior changes:
        Created variables, demon_cm, yWave_cm, xWave_cm, and two_pi to reduce the total amount of calculations neccessary.
        Moved the assignment for each index of amplitudeOut outside of the inner most for-loop and placed it after the 3rd for-loop. 

        New Changes:
        Removed denom_cm variable and eliminated using ints height and width because they are both equal to N. 
        Created local variables realOut_sum and imagOut_sum, so that the program would not have to reference 
        the same indecies of realOut and imagOut repeatedly inside the inner most for-loop.
*/
void fft_elim_mem_refs(){
    double realOut[N][N] = {0}, imagOut[N][N] = {0}, inputData[N][N] = {0}, amplitudeOut[N][N] = {0};
    //int height = N,  width = N; /* Unneccessary use of variables, no need for height and width since their values never change*/
    int yWave, xWave, ySpace, xSpace, i, j;
   //initialize the array inputData
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            #ifndef TEST
            inputData[i][j] = rand();
            #else
            inputData[i][j] = 1.0;
            #endif

    /* Saving the program from having to repeat calculations */
    double two_pi = 2 * PI;

    // Two outer loops iterate on output data.
    for (yWave = 0; yWave < N; yWave++) {
        /* Calculations replaced from the inner three loops */
        double yWave_cm = 1.0 * yWave / N; 

        for (xWave = 0; xWave < N; xWave++) {
            /* Accumulating result in local variables realOut_sum and imagOut_sum */
            double realOut_sum = 0, imagOut_sum = 0; 
            
            /* Calculations replaced from the inner two loops */
            double xWave_cm = 1.0 * xWave / N; 

            // Two inner loops iterate on input data.
            for (ySpace = 0; ySpace < N; ySpace++) {
                for (xSpace = 0; xSpace < N; xSpace++) {
                    // Compute real, imag, and ampltude.
                    double input_data = inputData[ySpace][xSpace]; /* Not on write up:  Created to save processor from having ot access this twice. */
                    realOut_sum += (input_data * cos(two_pi * ((1.0 * xWave_cm * xSpace) + (1.0 * yWave_cm * ySpace)))) / N; /* Uneccessary use of sqrt() function since the sqrt of N^2 will always be N. */
                    imagOut_sum -= (input_data * sin(two_pi * ((1.0 * xWave_cm * xSpace) + (1.0 * yWave_cm * ySpace)))) / N; /* Uneccessary use of sqrt() function since the sqrt of N^2 will always be N. */
                }
            }

            /* Calculating final value */
            realOut[yWave][xWave] = realOut_sum; /* Assigning acculumated value */
            imagOut[yWave][xWave] = imagOut_sum; /* Assigning acculumated value */
            amplitudeOut[yWave][xWave] = sqrt(realOut_sum * realOut_sum + imagOut_sum * imagOut_sum); /* Does not need to be in the inner for loop if called after */
        }
    }

    /* Checking to see if code matches */
    #ifdef TEST
    check_results(realOut, imagOut, amplitudeOut);
    #endif
}

/* FFT function – 4x1 loop unrolling. 
    Summary of changes:
        New changes:
            Implemented 4x1 loop unrolling. 
        
        Prior changes:
        Created variables, yWave_cm, xWave_cm, and two_pi to reduce the total amount of calculations neccessary. 
        Eliminated using ints height and width because they are both equal to N.
        Moved the assignment for each index of amplitudeOut outside of the inner most for-loop and placed it after the 3rd for-loop. 
        Created local variables realOut_sum and imagOut_sum, so that the program would not have to reference 
        the same indecies of realOut and imagOut repeatedly inside the inner most for-loop.
*/
void fft_loop_unrolling_4_1(){
    double realOut[N][N] = {0}, imagOut[N][N] = {0}, inputData[N][N] = {0}, amplitudeOut[N][N] = {0};
    int yWave, xWave, ySpace, xSpace, i, j;
   //initialize the array inputData
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            #ifndef TEST
            inputData[i][j] = rand();
            #else
            inputData[i][j] = 1.0;
            #endif
    /* Saving the program from having to repeat calculations */
    double two_pi = 2 * PI;
    // Two outer loops iterate on output data.
    for (yWave = 0; yWave < N; yWave++) {
        /* Calculations replaced from the inner three loops */
        double yWave_cm = 1.0 * yWave / N; 

        for (xWave = 0; xWave < N; xWave++) {
            /* Accumulating result in local variables realOut_sum and imagOut_sum */
            double realOut_sum = 0, imagOut_sum = 0; 
            
            /* Calculations replaced from the inner two loops */
            double xWave_cm = 1.0 * xWave / N; 
            
            /* Loop unrolling vars */
            long limit = N - 3; 
            double x = IDENT, y = IDENT;

            // Two inner loops iterate on input data.
            for (ySpace = 0; ySpace < N; ySpace++) {
                for (xSpace = 0; xSpace < limit; xSpace+=4) {
                    /* Not on write up:  Created to save processor from having ot access this twice. */
                    double input_data0 = inputData[ySpace][xSpace]; 
                    double input_data1 = inputData[ySpace][xSpace + 1]; 
                    double input_data2 = inputData[ySpace][xSpace + 2]; 
                    double input_data3 = inputData[ySpace][xSpace + 3]; 

                    // Compute real, imag, and ampltude.
                    x += ((input_data0 * cos(two_pi * ((1.0 * xWave_cm * xSpace) + (1.0 * yWave_cm * ySpace)))) / N) 
                        OP ((input_data1 * cos(two_pi * ((1.0 * xWave_cm * (xSpace + 1)) + (1.0 * yWave_cm * ySpace)))) / N)
                        OP ((input_data2 * cos(two_pi * ((1.0 * xWave_cm * (xSpace + 2)) + (1.0 * yWave_cm * ySpace)))) / N)
                        OP ((input_data3 * cos(two_pi * ((1.0 * xWave_cm * (xSpace + 3)) + (1.0 * yWave_cm * ySpace)))) / N); 
                    y -= ((input_data0 * sin(two_pi * ((1.0 * xWave_cm * xSpace) + (1.0 * yWave_cm * ySpace)))) / N) 
                        OP ((input_data1 * sin(two_pi * ((1.0 * xWave_cm * (xSpace + 1)) + (1.0 * yWave_cm * ySpace)))) / N)
                        OP ((input_data2 * sin(two_pi * ((1.0 * xWave_cm * (xSpace + 2)) + (1.0 * yWave_cm * ySpace)))) / N)
                        OP ((input_data3 * sin(two_pi * ((1.0 * xWave_cm * (xSpace + 3)) + (1.0 * yWave_cm * ySpace)))) / N); 
                }
                for (; xSpace < N; xSpace++) {
                    /* Not on write up:  Created to save processor from having ot access this twice. */
                    double input_data0 = inputData[ySpace][xSpace]; 

                    x += (input_data0 * cos(two_pi * ((1.0 * xWave_cm * xSpace) + (1.0 * yWave_cm * ySpace)))) / N; /* Uneccessary use of sqrt() function since the sqrt of N^2 will always be N. */
                    y -= (input_data0 * sin(two_pi * ((1.0 * xWave_cm * xSpace) + (1.0 * yWave_cm * ySpace)))) / N; /* Uneccessary use of sqrt() function since the sqrt of N^2 will always be N. */
                }
            }

            /* Calculating final value */
            realOut[yWave][xWave] = x; /* Assigning acculumated value */
            imagOut[yWave][xWave] = y; /* Assigning acculumated value */
            amplitudeOut[yWave][xWave] = sqrt(x * x + y * y); /* Does not need to be in the inner for loop if called after */
        }
    }

    /* Checking to see if code matches */
    #ifdef TEST
    check_results(realOut, imagOut, amplitudeOut);
    #endif
}

/* FFT function – 8x1 loop unrolling. 
    Summary of changes:
        New changes:
            Implemented 8x1 loop unrolling. 
        
        Prior changes:
        Created variables, yWave_cm, xWave_cm, and two_pi to reduce the total amount of calculations neccessary. 
        Eliminated using ints height and width because they are both equal to N.
        Moved the assignment for each index of amplitudeOut outside of the inner most for-loop and placed it after the 3rd for-loop. 
        Created local variables realOut_sum and imagOut_sum, so that the program would not have to reference 
        the same indecies of realOut and imagOut repeatedly inside the inner most for-loop.
*/
void fft_loop_unrolling_8_1(){
    double realOut[N][N] = {0}, imagOut[N][N] = {0}, inputData[N][N] = {0}, amplitudeOut[N][N] = {0};
    int yWave, xWave, ySpace, xSpace, i, j;
   //initialize the array inputData
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            #ifndef TEST
            inputData[i][j] = rand();
            #else
            inputData[i][j] = 1.0;
            #endif
    /* Saving the program from having to repeat calculations */
    double two_pi = 2 * PI;
    // Two outer loops iterate on output data.
    for (yWave = 0; yWave < N; yWave++) {
        /* Calculations replaced from the inner three loops */
        double yWave_cm = 1.0 * yWave / N; 

        for (xWave = 0; xWave < N; xWave++) {
            /* Accumulating result in local variables realOut_sum and imagOut_sum */
            double realOut_sum = 0, imagOut_sum = 0; 
            
            /* Calculations replaced from the inner two loops */
            double xWave_cm = 1.0 * xWave / N; 

            /* Loop unrolling vars */
            long limit = N - 7; 
            double x = IDENT, y = IDENT;

            // Two inner loops iterate on input data.
            for (ySpace = 0; ySpace < N; ySpace++) {
                for (xSpace = 0; xSpace < limit; xSpace+=8) {
                    /* Not on write up:  Created to save processor from having ot access this twice. */
                    double input_data0 = inputData[ySpace][xSpace]; 
                    double input_data1 = inputData[ySpace][xSpace + 1]; 
                    double input_data2 = inputData[ySpace][xSpace + 2]; 
                    double input_data3 = inputData[ySpace][xSpace + 3]; 
                    double input_data4 = inputData[ySpace][xSpace + 4]; 
                    double input_data5 = inputData[ySpace][xSpace + 5]; 
                    double input_data6 = inputData[ySpace][xSpace + 6]; 
                    double input_data7 = inputData[ySpace][xSpace + 7]; 
                    
                    // Compute real, imag, and ampltude.
                    x += ((input_data0 * cos(two_pi * ((1.0 * xWave_cm * xSpace) + (1.0 * yWave_cm * ySpace)))) / N) 
                        OP ((input_data1 * cos(two_pi * ((1.0 * xWave_cm * (xSpace + 1)) + (1.0 * yWave_cm * ySpace)))) / N)
                        OP ((input_data2 * cos(two_pi * ((1.0 * xWave_cm * (xSpace + 2)) + (1.0 * yWave_cm * ySpace)))) / N)
                        OP ((input_data3 * cos(two_pi * ((1.0 * xWave_cm * (xSpace + 3)) + (1.0 * yWave_cm * ySpace)))) / N)
                        OP ((input_data4 * cos(two_pi * ((1.0 * xWave_cm * (xSpace + 4)) + (1.0 * yWave_cm * ySpace)))) / N)
                        OP ((input_data5 * cos(two_pi * ((1.0 * xWave_cm * (xSpace + 5)) + (1.0 * yWave_cm * ySpace)))) / N)
                        OP ((input_data6 * cos(two_pi * ((1.0 * xWave_cm * (xSpace + 6)) + (1.0 * yWave_cm * ySpace)))) / N)
                        OP ((input_data7 * cos(two_pi * ((1.0 * xWave_cm * (xSpace + 7)) + (1.0 * yWave_cm * ySpace)))) / N); 
                    y -= ((input_data0 * sin(two_pi * ((1.0 * xWave_cm * xSpace) + (1.0 * yWave_cm * ySpace)))) / N) 
                        OP ((input_data1 * sin(two_pi * ((1.0 * xWave_cm * (xSpace + 1)) + (1.0 * yWave_cm * ySpace)))) / N)
                        OP ((input_data2 * sin(two_pi * ((1.0 * xWave_cm * (xSpace + 2)) + (1.0 * yWave_cm * ySpace)))) / N)
                        OP ((input_data3 * sin(two_pi * ((1.0 * xWave_cm * (xSpace + 3)) + (1.0 * yWave_cm * ySpace)))) / N)
                        OP ((input_data4 * sin(two_pi * ((1.0 * xWave_cm * (xSpace + 4)) + (1.0 * yWave_cm * ySpace)))) / N)
                        OP ((input_data5 * sin(two_pi * ((1.0 * xWave_cm * (xSpace + 5)) + (1.0 * yWave_cm * ySpace)))) / N)
                        OP ((input_data6 * sin(two_pi * ((1.0 * xWave_cm * (xSpace + 6)) + (1.0 * yWave_cm * ySpace)))) / N)
                        OP ((input_data7 * sin(two_pi * ((1.0 * xWave_cm * (xSpace + 7)) + (1.0 * yWave_cm * ySpace)))) / N);
                }
                for (; xSpace < N; xSpace++) {
                    /* Not on write up:  Created to save processor from having ot access this twice. */
                    double input_data0 = inputData[ySpace][xSpace]; 

                    x += (input_data0 * cos(two_pi * ((1.0 * xWave_cm * xSpace) + (1.0 * yWave_cm * ySpace)))) / N; /* Uneccessary use of sqrt() function since the sqrt of N^2 will always be N. */
                    y -= (input_data0 * sin(two_pi * ((1.0 * xWave_cm * xSpace) + (1.0 * yWave_cm * ySpace)))) / N; /* Uneccessary use of sqrt() function since the sqrt of N^2 will always be N. */
                }
            }

            /* Calculating final value */
            realOut[yWave][xWave] = x; /* Assigning acculumated value */
            imagOut[yWave][xWave] = y; /* Assigning acculumated value */
            amplitudeOut[yWave][xWave] = sqrt(x * x + y * y); /* Does not need to be in the inner for loop if called after */
        }
    }

    /* Checking to see if code matches */
    #ifdef TEST
    check_results(realOut, imagOut, amplitudeOut);
    #endif
}

/* FFT function – 2x1 loop unrolling. 
    Summary of changes:
        New changes:
            Implemented 2x1 loop unrolling. 
        
        Prior changes:
        Created variables, yWave_cm, xWave_cm, and two_pi to reduce the total amount of calculations neccessary. 
        Eliminated using ints height and width because they are both equal to N.
        Moved the assignment for each index of amplitudeOut outside of the inner most for-loop and placed it after the 3rd for-loop. 
        Created local variables realOut_sum and imagOut_sum, so that the program would not have to reference 
        the same indecies of realOut and imagOut repeatedly inside the inner most for-loop.
*/
void fft_loop_unrolling_2_1(){
    double realOut[N][N] = {0}, imagOut[N][N] = {0}, inputData[N][N] = {0}, amplitudeOut[N][N] = {0};
    int yWave, xWave, ySpace, xSpace, i, j;
   //initialize the array inputData
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            #ifndef TEST
            inputData[i][j] = rand();
            #else
            inputData[i][j] = 1.0;
            #endif
    /* Saving the program from having to repeat calculations */
    double two_pi = 2 * PI;
    // Two outer loops iterate on output data.
    for (yWave = 0; yWave < N; yWave++) {
        /* Calculations replaced from the inner three loops */
        double yWave_cm = 1.0 * yWave / N; 

        for (xWave = 0; xWave < N; xWave++) {
            /* Accumulating result in local variables realOut_sum and imagOut_sum */
            double realOut_sum = 0, imagOut_sum = 0; 
            
            /* Calculations replaced from the inner two loops */
            double xWave_cm = 1.0 * xWave / N; 

            /* Loop unrolling vars */
            long limit = N - 1; 
            double x = IDENT, y = IDENT;

            // Two inner loops iterate on input data.
            for (ySpace = 0; ySpace < N; ySpace++) {
                for (xSpace = 0; xSpace < limit; xSpace+=2) {
                    /* Not on write up:  Created to save processor from having ot access this twice. */
                    double input_data0 = inputData[ySpace][xSpace]; 
                    double input_data1 = inputData[ySpace][xSpace + 1]; 
                    // Compute real, imag, and ampltude.
                    x += ((input_data0 * cos(two_pi * ((1.0 * xWave_cm * xSpace) + (1.0 * yWave_cm * ySpace)))) / N) 
                        OP ((input_data1 * cos(two_pi * ((1.0 * xWave_cm * (xSpace + 1)) + (1.0 * yWave_cm * ySpace)))) / N); 
                    y -= ((input_data0 * sin(two_pi * ((1.0 * xWave_cm * xSpace) + (1.0 * yWave_cm * ySpace)))) / N) 
                        OP ((input_data1 * sin(two_pi * ((1.0 * xWave_cm * (xSpace + 1)) + (1.0 * yWave_cm * ySpace)))) / N); 
                }
                for (; xSpace < N; xSpace++) {
                    /* Not on write up:  Created to save processor from having ot access this twice. */
                    double input_data0 = inputData[ySpace][xSpace]; 

                    x += (input_data0 * cos(two_pi * ((1.0 * xWave_cm * xSpace) + (1.0 * yWave_cm * ySpace)))) / N; /* Uneccessary use of sqrt() function since the sqrt of N^2 will always be N. */
                    y -= (input_data0 * sin(two_pi * ((1.0 * xWave_cm * xSpace) + (1.0 * yWave_cm * ySpace)))) / N; /* Uneccessary use of sqrt() function since the sqrt of N^2 will always be N. */
                }
            }

            /* Calculating final value */
            realOut[yWave][xWave] = x; /* Assigning acculumated value */
            imagOut[yWave][xWave] = y; /* Assigning acculumated value */
            amplitudeOut[yWave][xWave] = sqrt(x * x + y * y); /* Does not need to be in the inner for loop if called after */
        }
    }

    /* Checking to see if code matches */
    #ifdef TEST
    check_results(realOut, imagOut, amplitudeOut);
    #endif
}

/* FFT function – 2x2 loop unrolling. 
    Summary of changes:
        New changes:
            Implemented 2x2 loop unrolling. 
        
        Prior changes:
        Created variables, yWave_cm, xWave_cm, and two_pi to reduce the total amount of calculations neccessary. 
        Eliminated using ints height and width because they are both equal to N.
        Moved the assignment for each index of amplitudeOut outside of the inner most for-loop and placed it after the 3rd for-loop. 
        Created local variables realOut_sum and imagOut_sum, so that the program would not have to reference 
        the same indecies of realOut and imagOut repeatedly inside the inner most for-loop.
*/
void fft_loop_unrolling_2_2(){
    double realOut[N][N] = {0}, imagOut[N][N] = {0}, inputData[N][N] = {0}, amplitudeOut[N][N] = {0};
    //int height = N,  width = N; /* Unneccessary use of variables, no need for height and width since their values never change*/
    int yWave, xWave, ySpace, xSpace, i, j;
   //initialize the array inputData
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            #ifndef TEST
            inputData[i][j] = rand();
            #else
            inputData[i][j] = 1.0;
            #endif
    /* Saving the program from having to repeat calculations */
    double two_pi = 2 * PI;
    // Two outer loops iterate on output data.
    for (yWave = 0; yWave < N; yWave++) {
        /* Calculations replaced from the inner three loops */
        double yWave_cm = 1.0 * yWave / N; 

        for (xWave = 0; xWave < N; xWave++) {
            /* Accumulating result in local variables realOut_sum and imagOut_sum */
            double realOut_sum = 0, imagOut_sum = 0; 
            
            /* Calculations replaced from the inner two loops */
            double xWave_cm = 1.0 * xWave / N; 
            


            /* Loop unrolling vars */
            long limit = N - 1; 
            double x0 = IDENT, x1 = IDENT, y0 = IDENT, y1 = IDENT; /* Copied from L10 Slides */

            // Two inner loops iterate on input data.
            for (ySpace = 0; ySpace < N; ySpace++) {
                for (xSpace = 0; xSpace < limit; xSpace+=2) {
                    /* Not on write up:  Created to save processor from having ot access this twice. */
                    double input_data0 = inputData[ySpace][xSpace]; 
                    double input_data1 = inputData[ySpace][xSpace + 1]; 
                    // Compute real, imag, and ampltude.
                    x0 += ((input_data0 * cos(two_pi * ((1.0 * xWave_cm * xSpace) + (1.0 * yWave_cm * ySpace)))) / N);
                    x1 += ((input_data1 * cos(two_pi * ((1.0 * xWave_cm * (xSpace + 1)) + (1.0 * yWave_cm * ySpace)))) / N); 
                    y0 -= ((input_data0 * sin(two_pi * ((1.0 * xWave_cm * xSpace) + (1.0 * yWave_cm * ySpace)))) / N);
                    y1 -= ((input_data1 * sin(two_pi * ((1.0 * xWave_cm * (xSpace + 1)) + (1.0 * yWave_cm * ySpace)))) / N); 
                }
                for (; xSpace < N; xSpace++) {
                    double input_data0 = inputData[ySpace][xSpace]; 
                    
                    x0 += (input_data0 * cos(two_pi * ((1.0 * xWave_cm * xSpace) + (1.0 * yWave_cm * ySpace)))) / N; /* Uneccessary use of sqrt() function since the sqrt of N^2 will always be N. */
                    y0 -= (input_data0 * sin(two_pi * ((1.0 * xWave_cm * xSpace) + (1.0 * yWave_cm * ySpace)))) / N; /* Uneccessary use of sqrt() function since the sqrt of N^2 will always be N. */
                }
            }

            /* Calculating final value */
            realOut[yWave][xWave] = x0 OP x1; /* Assigning acculumated value */
            imagOut[yWave][xWave] = y0 OP y1; /* Assigning acculumated value */
            amplitudeOut[yWave][xWave] = sqrt(realOut[yWave][xWave] * realOut[yWave][xWave] + imagOut[yWave][xWave] * imagOut[yWave][xWave]); /* Does not need to be in the inner for loop if called after */
        }
    }

    /* Checking to see if code matches */
    #ifdef TEST
    check_results(realOut, imagOut, amplitudeOut);
    #endif
}

/* Floating point equal function - Helper function that accounts for rounding error */
int floating_point_equal(double a, double b){
    return fabs(a-b) < 0.00000000001; /*  If I add anymore precision I get roudning errors */
}

/* Check results function – checks to see if code of 2D arrays match. */
void check_results(double realOut[N][N], double imagOut[N][N], double amplitudeOut[N][N]) {
    for (int x = 0; x < N; x++) {
        for (int y = 0; y < N; y++) {
            if ((floating_point_equal(realOut_ref[x][y],realOut[x][y])) == 0) {
                printf("Failed\n\nIndex (%d, %d): %lf did not match %lf\t", x, y, realOut[x][y], realOut_ref[x][y]);
                printf("- Failed REAL test.\n");
                return;
            } else if ((floating_point_equal(imagOut_ref[x][y],imagOut[x][y])) == 0) {
                printf("Failed\n\nIndex (%d, %d): %lf did not match %lf\t", x, y, imagOut[x][y], imagOut_ref[x][y]);
                printf("- Failed IMAG test.\n");
                return;
            } else if ((floating_point_equal(amplitudeOut_ref[x][y],amplitudeOut[x][y])) == 0) {
                printf("Failed\n\nIndex (%d, %d): %lf did not match %lf\t", x, y, amplitudeOut[x][y], amplitudeOut_ref[x][y]);
                printf("- Failed AMP test.\n");
                return;
            }
        }
    }

    printf("passed\t\t");
}