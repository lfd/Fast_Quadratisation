###
Create and run the docker file:
```
docker build --no-cache -t ext .

docker run -it ext
```
Run the experiments:

```
source activate quark_install

cd reduc/

python Test.py
```

Create Plots:
```
Rscript plots.r
```

Export Plots:
```
docker cp <container_name>:/repro/ <host_file_system>
```