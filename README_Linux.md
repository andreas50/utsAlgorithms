### Compile demo

```
gcc -Wall ema.c sma.c rolling.c test.c -o test -lm
./test
```

### Generate dynamically linked shared object library

```
gcc -Wall -fPIC -shared sma.c ema.c rolling.c -o libUTSOperators.so
```

### Compile demo via shared library

Linking

```
gcc -Wall -L./ test.c -o test -lUTSOperators -lm
```

Run demo (shared library made available via `LD_LIBRARY_PATH`)

```
export LD_LIBRARY_PATH=`pwd`:$LD_LIBRARY_PATH
echo $LD_LIBRARY_PATH
./test
unset LD_LIBRARY_PATH
```

Run demo (shared library made available by copying to standard library location)

```
sudo cp libUTSOperators.so /usr/lib
sudo chmod 755 /usr/lib/libUTSOperators.so  # make library available to all users
ldd test                                    # check that liboperators.so can be found for execution
./test
sudo rm /usr/lib/libUTSOperators.so
```
