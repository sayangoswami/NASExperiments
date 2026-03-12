# Selective Sequencing Experiments


## Setting up

1. Install conda, if not already installed, [from here](https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html#installing-conda).
2. Create a conda environment and activate.
	```bash
	conda env create -f env.yml
	conda activate selectiveseq
	```
3. Install Metagraph prerequisites [from here](https://metagraph.ethz.ch/static/docs/installation.html#prerequisites).
4. Clone repo and setup.
	```bash
	git clone <experiments repo>
	cd NASExperiments
	chmod +x setup.sh
	./setup.sh
	```
5.  By default, `ont-pybasecall-client-lib=7.11.2` and `ont-pyguppy-client-lib=6.5.7` are installed in the environment. To check if they are compatible with the basecaller vesion, run `/opt/ont/dorado/bin/dorado_basecall_server --version` and, if required, download the appropriate client library from [here if you have Dorado 7.3.0 onwards](https://pypi.org/project/ont-pybasecall-client-lib/) or [here if you have Guppy or Dorado upto version 7.2.x](https://pypi.org/project/ont-pyguppy-client-lib/).


## Running experiments on basecalled reads

In the `experiments` folder, so the following:

1. Create a script (say `my_tasks.sh`) as follows. Use the variables `OUTDIR` for output files and `TMPDIR` for temporary files. This script will be imported by another script (`task_runner.sh`) and these variables would be defined in it. The functions should be called `fn` (n = 1, 2, 3, ...) and define `ALL_TASKS="1-n"`. Collinearity, Metagraph, Spumoni, and Rawhash executable paths are all exported in `task_runner.sh`.
	
```bash
ALL_TASKS="1-3"

f1() {
	# task 1
	Collinearity ...
}

f2() {
	# task 2
	metagraph ...
}

f3() {
	# task 3
	spumoni ...
}
```
	
2. Run the tasks using task runner as follows:

```bash
# Run all tasks from 'my_tasks.sh'
./task_runner.sh -t my_tasks.sh

# Can also omit the .sh extension from task file
./task-runner.sh -t my_tasks

# Run tasks 1-2 from the 'my_tasks.sh' file
./task_runner.sh -t my_tasks 1-2

# Run ./task_runner.sh -h for other formats to specify task numbers
```

### Logs

- The stdout and stderr logs are written in a file `{TASK_NAME}_MMMDD_HHMMSS.log`. 
- The output directory variable `OUTDIR` points to `NASExperiments/out`
- `TMPDIR` points to `NASExperiments/tmp`
## Running experiments on raw signals using Minknow Simulator

Run the server -

```bash
mksimserver --certs /scratch/NASExperiments/code/MinknoApiSimulator/certs \
--input /data/SimulatedDatasets/Zymo/signals/Sigs0_450.blow5 \
--input /data/SimulatedDatasets/Zymo/signals/Sigs1_450.blow5
```

On a different shell, run Readfish (assuming the `.toml` file is validated using `readfish validate`)-

```bash
export MINKNOW_API_USE_LOCAL_TOKEN="no"
export MINKNOW_SIMULATOR="true"
export MINKNOW_TRUSTED_CA=/scratch/NASExperiments/code/MinknoApiSimulator/certs/server.pem
export MINKNOW_API_CLIENT_CERTIFICATE_CHAIN=/scratch/NASExperiments/code/MinknoApiSimulator/certs/client.pem
export MINKNOW_API_CLIENT_KEY=/scratch/NASExperiments/code/MinknoApiSimulator/certs/client.key
 
readfish targets --wait-for-ready 5 \
--toml /scratch/NASExperiments/configs/rf_mm_zymo.toml \
--port 50051 --device MN12345 \
--log-file /scratch/NASExperiments/logs/readfish_test.log \
--experiment-name 'test.log'
```

Or run the `NASExperiments/scripts/simulate_run.sh` script after editing the following lines:

```bash
INPUT_SIGNAL=(
	...
)
CONFIG_TOML=...
```

### Logs

Readfish creates a `test_run_readfish.tsv` file in the directory where its invoked. It contains the decision taken for each read. Two additional files are created in `NASExperiments/logs` - `readfish_MMMDD_hhmmss.log` and `server_MMMDD_hhmmss.log` which contains the stdout and stderr logs for `readfish` and `mksimserver` respectively.
