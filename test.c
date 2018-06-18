// Copyright 2012-2017 by Andreas Eckner
// License GPL-2 | GPL-3

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "ema.h"
#include "sma.h"
#include "rolling.h"


// Print nicely formatted observation times and values for an unevenly spaced time series
void print_uts(double values[], double times[], int n)
{
  // values     ... array of observation values
  // times      ... array of observation times
  // n          ... number of observations, i.e. length of 'values' and 'times'  
  
  // Print header
  printf("-------------\n");
  printf("Time    Value\n");
  printf("-------------\n");

  for (int i=0; i < n; i++)
    printf("%.1f %9.2f\n", times[i], values[i]); 
}


// Demo of functionality
int main()
{
  // Define sample time series
  double values[] = {0, 2, 4, 6, 8, 10};
  double times[] = {0, 1, 1.2, 2.3, 2.9, 5};
  int n = sizeof(values) / sizeof(double);
  double out[n];
  printf("Input time series X\n");
  print_uts(values, times, n);
  
  // Parameters for the demo
  double width_before = 2.5, width_after = 1, width_before_long = 1000;
  double tau = 1.5, tau_long = 1000;

  /*
    Basic Rolling Time Series Operators
  */
  printf("\n\n##### Basic Rolling Time Series Operators #####\n");

  // rolling numer of observations
  rolling_num_obs(values, times, &n, out, &width_before, &width_after);
  printf("\nrolling_num_obs(X, %.1f, %.1f)\n", width_before, width_after);
  print_uts(out, times, n);

  // rolling sum
  rolling_sum(values, times, &n, out, &width_before, &width_after);
  printf("\nrolling_sum(X, %.1f, %.1f)\n", width_before, width_after);
  print_uts(out, times, n);

  // rolling product
  rolling_product(values, times, &n, out, &width_before, &width_after);
  printf("\nrolling_product(X, %.1f, %.1f)\n", width_before, width_after);
  print_uts(out, times, n);

  // rolling average
  rolling_mean(values, times, &n, out, &width_before, &width_after);
  printf("\nrolling_mean(X, %.1f, %.1f)\n", width_before, width_after);
  print_uts(out, times, n);

  // rolling median
  rolling_median(values, times, &n, out, &width_before, &width_after);
  printf("\nrolling_median(X, %.1f, %.1f)\n", width_before, width_after);
  print_uts(out, times, n);
  
  // rolling maximum
  rolling_max(values, times, &n, out, &width_before, &width_after);
  printf("\nrolling_max(X, %.1f, %.1f)\n", width_before, width_after);
  print_uts(out, times, n);
  
  // rolling minimum
  rolling_min(values, times, &n, out, &width_before, &width_after);
  printf("\nrolling_min(X, %.1f, %.1f)\n", width_before, width_after);
  print_uts(out, times, n);
  
  // rolling standard deviation
  rolling_sd(values, times, &n, out, &width_before, &width_after);
  printf("\nrolling_sd(X, %.1f, %.1f)\n", width_before, width_after);
  print_uts(out, times, n);
  
  // rolling variance
  rolling_var(values, times, &n, out, &width_before, &width_after);
  printf("\nrolling_var(X, %.1f, %.1f)\n", width_before, width_after);
  print_uts(out, times, n);


  /*
    Simple Moving Averages (SMAs)
  */
  printf("\n\n##### Simple Moving Averages (SMAs) #####\n");

  // SMA with last-point interpolation
  sma_last(values, times, &n, out, &width_before, &width_after);
  printf("\nSMA_last(X, %.1f, %.1f)\n", width_before, width_after);
  print_uts(out, times, n);
  
  // SMA with next-point interpolation
  sma_next(values, times, &n, out, &width_before, &width_after);
  printf("\nSMA_next(X, %.1f, %.1f)\n", width_before, width_after);
  print_uts(out, times, n);

  // SMA with linear interpolation
  sma_linear(values, times, &n, out, &width_before, &width_after);
  printf("\nSMA_linear(X, %.1f, %.1f)\n", width_before, width_after);
  print_uts(out, times, n);
  
  // SMA with slow time decay
  sma_linear(values, times, &n, out, &width_before_long, &width_after);
  printf("\nSMA_linear(X, %.1f, %.1f) ... a SMA with a long rolling time window produces nearly constant output\n", width_before_long, width_after);
  print_uts(out, times, n);
  

  /*
    Explonential Moving Averages (EMAs)
  */
  printf("\n\n##### Exponential Moving Averages (EMAs) #####\n");

  // EMA with last-point interpolation
  ema_last(values, times, &n, out, &tau);
  printf("\nEMA_last(X, %.1f)\n", tau);
  print_uts(out, times, n);

  // EMA with next-point interpolation
  ema_next(values, times, &n, out, &tau);
  printf("\nEMA_next(X, %.1f)\n", tau);
  print_uts(out, times, n);

  // EMA with linear interpolation
  ema_linear(values, times, &n, out, &tau);
  printf("\nEMA_linear(X, %.1f)\n", tau);
  print_uts(out, times, n);
  
  // EMA with slow time decay
  ema_linear(values, times, &n, out, &tau_long);
  printf("\nEMA_linear(X, %.1f) ... an EMA with a slow time decay produces nearly constant output\n", tau_long);
  print_uts(out, times, n);

  // Wait for key pressed before exiting
  printf("\nPress <ENTER> to exit the program.\n");
  getchar();
  return 0;
}
