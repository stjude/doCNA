# doCNA

Tool for analysis of copy number profile in WGS.

## Installing doCNA

```shell
git clone https://github.com/stjude/doCNA.git

# setup a virtual environment
python -m venv docna_install
source ./docna_install/bin/activate

# install into virtual environment
pip3 install ./doCNA
```

# download the appropriate reference files from Zenodo (hg19 and hg38 available)
[Zenodo Refs](https://sandbox.zenodo.org/record/1181478#.ZDhvH7qZNPY)
```
# hg19 example
wget https://sandbox.zenodo.org/record/1181478/files/hg19_cytoBand.dat && echo 4036fc07b0b87ef28b46b1229141abb1 hg19_cytoBand.dat | md5sum --check
...
hg19_cytoBand.dat: OK

wget https://sandbox.zenodo.org/record/1181478/files/hg19_SuperGood.dat.gz && echo 8e3058f18c502a91b466d779a59a35f0 hg19_SuperGood.dat.gz | md5sum --check
...
hg19_SuperGood.dat.gz: OK

# hg38 example
wget https://sandbox.zenodo.org/record/1181478/files/hg38_cytoBand.dat && echo 5c957c934461320fdf6211df3d68bdd3 hg38_cytoBand.dat | md5sum --check
...
hg38_cytoBand.dat: OK


wget https://sandbox.zenodo.org/record/1181478/files/hg38_SuperGood.dat.gz && echo 1c32772977772ee1bab65f4be3a7acf2 hg38_SuperGood.dat.gz | md5sum --check
...
hg38_SuperGood.dat.gz: OK
```

## Workflow for running doCNA

1. Get a local copy of the config:
```shell
docna getconfig
```

2. Edit the config to point to real SuperGood and CytoBand files. We don't recommend changing additional parameters when starting out.

3. Run doCNA to generate output files:
```shell
# replace sample with whatever your samplename is

docna analyze -i Sample.txt -c doCNA.ini -s Sample
... Lots of log info! ...
12:18:19 doCNA.WGS: INFO: Ready to report!
All done
```

4. Launch the viewer and visualize output:
```shell
# if you are running this on your local machine, such as a desktop or laptop:
docna viewer

# if you are running this on a remote machine, like a cluster:
docna viewer --remote
**********
Access dashboard in browser via: http://10.220.16.129:39738
**********
INFO:     Started server process [11716]
INFO:     Waiting for application startup.
INFO:     Application startup complete.
INFO:     Uvicorn running on http://10.220.16.129:39738 (Press CTRL+C to quit)
```
5. Load files:
![example](./examples/docna.gif)
