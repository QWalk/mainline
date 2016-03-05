import os
import shutil as shu

def fix_traces():
  """Ensure all trace files have correct naming."""
  files = os.listdir("./")
  traces = [f for f in files if "trace" in f]
  bad_traces = [f for f in traces if "tmoves" not in f]
  if len(bad_traces)==0:
    print "Fix trace: Found no bad traces."
    return None
  dmcs = [f for f in files if f[-3:] == "dmc"]
  base = dmcs[0].split('_')[0]
  fix_traces = [ trace.replace("qw",base) for trace in bad_traces ]
  for i in range(len(bad_traces)):
    shu.copy(bad_traces[i],fix_traces[i])
    print "Fix trace: copying %s to %s"%(
        bad_traces[i],fix_traces[i]
      )
