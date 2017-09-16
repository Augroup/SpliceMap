# SpliceMap
A fork of Standford's SpliceMap -- Copyright (c) 2010 Kin Fai Au and John Chong Mu.

## How to run from the docker

#### Testing the command:
The docker image comes with a small example data.  You can access this and run the command on the example which will create an `output` folder and a `temp` folder in your current directory.

```bash
$ docker run -v $(pwd)/output:/Source/SpliceMap/example/output \
             -v $(pwd)/temp:/Source/SpliceMap/example/temp \
             -w /Source/SpliceMap/example \
             vacation/SpliceMap runSpliceMap run.cfg
```

#### Running the command on your data:

Construct a `run.cfg` pointing to required files within your working directory such as this one:

https://github.com/jason-weirather/SpliceMap/tree/master/example

And from within your working directory you can call splicemap on your `run.cfg` file like so:

```bash
$ docker run -v $(pwd):/home vacation/SpliceMap runSpliceMap run.cfg
```
