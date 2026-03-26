## Reproduction package for [It’s Quick to be Square: Fast Quadratisation for Quantum Toolchains](https://doi.org/10.1145/3800943)

### Create and run the docker file:
```
docker build --no-cache -t ext .

docker run -it ext
```

### Extract experiment data:
```
make extract
```

### Run the experiments:

```
source activate quark_install

make data
```


### Create Plots:
```
make plots
```

### Export Plots:
```
docker cp <container_name>:/repro/ <host_file_system>
```