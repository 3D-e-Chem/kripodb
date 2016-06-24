# Kripo pipelines as CWL workflow

For more info on Common Workflow Language (CWL) see http://www.commonwl.org/

Requires CWL runner which can be installed with
```
pip install cwl-runner
```
(Note! CWL runner does not yet support Python3 see https://github.com/common-workflow-language/cwltool/issues/8)

## Compute distance matrix of Kripo fingerprints

The pipeline can be run on the tiny example dataset included in `/data` directory of this repo.

Prepare input files
```
cp ../../data/fragments.sqlite .
cwl-runner kripodb-fingerprints-export.cwl \
  --fingerprintdb ../../data/fingerprints.sqlite \
  --fingerprinttxt_name fingerprints.txt
```

Run workflow
```
cwl-runner --verbose kripo-fingerprints2matrix.cwl kripo-fingerprints2matrix.json
```

Check output using, pytables (http://www.pytables.org/) dumper
```
docker run \
  -v $PWD/dist.packedfrozen.h5:/data/dist.packedfrozen.h5:ro \
  rcsa/python2-hdf5 \
  ptdump -vd /data/dist.packedfrozen.h5 | less
```
The '/scores' node should contain a 1000x1000 matrix with some non-zero values.
