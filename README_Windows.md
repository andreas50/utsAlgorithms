### Compiler Setup

1.) Install MinGW

2.) Set up Windows path environment
```path = %PATH%;c:\mingw\bin```

3.) Go to the directory with the `C` source files



### Compile demo

gcc -std=c99 -Wall ema.c sma.c rolling.c test.c -o test -lm
test


##########################################################
# Generate DLL file and compile demo using this DLL file #
##########################################################

# Create DLL file
gcc -std=c99 -Wall -shared sma.c ema.c rolling.c -o UTSOperators.dll

# Compile demo against DLL file
gcc -std=c99 -L. -l UTSOperators test.c -o test
test
